package controllers

import utils.Utils._
import utils.BamUtils._
import utils.Constants._

import javax.inject.Inject
import java.io.{File, FileInputStream}
import java.nio.file.{Files, Path, Paths}

import play.api.mvc._
import play.api.libs.json.Json._
import play.api.db._
import play.api.Logger
import play.api.http.{HttpEntity, MimeTypes}

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.util.Random
import akka.stream.scaladsl.StreamConverters

import htsjdk.samtools.{SamReader, SamReaderFactory}
//import htsjdk.samtools.{SAMFileHeader, SAMFileWriter, SAMFileWriterFactory}
import sys.process._


class BamReader(val bam:String){
  val reader:SamReader = SamReaderFactory.makeDefault().open(new File(bam))
  def head() = {
    var i = 0
    val iterate = reader.iterator()
    while(iterate.hasNext && i < 10){
      println(iterate.next())
      i += 1
    }
  }
}


/**
  * Created by jdelafon on 06/03/16.
  */
class BamController @Inject()(db: Database) extends Controller {

  /*
   * Return the file names scorresponding to this key in the database
   */
  def getBamNames(key: String) : List[String] = {
    db.withConnection { conn =>
      val cursor = conn.createStatement
      val res = cursor.executeQuery("SELECT `key`,`filename` FROM bam WHERE `key`='"+key+"';")
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      Logger.info("File names: " + filenames.mkString)
      filenames.toList
    }
  }


  /*
   * Return the BAM text as output by "samtools view -hb <filename> <region>".
   * It is encoded in base64 because transfers have to be ascii (?).
   */
  def read(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamNames(key).head
    val bamPath: Path = Paths.get(BAM_PATH, filename)
    val command = s"samtools view -hb $bamPath $region" #| "base64"
    val res: String = command.!!
    // TODO: check the format of $region
    Ok(res).as(HTML)

    //import play.api.libs.streams.Streams
    //import akka.util.ByteString
    // import java.io.ByteArrayInputStream

    //val bam = bamPath.toFile
    //val stream: FileInputStream = new FileInputStream(bam)
    //val source = StreamConverters.fromInputStream(() => stream)
    //val contentLength = bam.length
    //Result(
    //  header = ResponseHeader(OK, Map(
    //    CONTENT_TYPE -> "application/octet-stream",
    //    CONTENT_LENGTH -> contentLength.toString,
    //    CONNECTION -> "keep-alive"
    //  )),
    //  body = HttpEntity.Streamed(source, Some(contentLength), Some("application/octet-stream"))
    //)
  }


  /*
   * Return a JSON summarizing the content of this BAM region (1 item per alignment)
   */
  def view(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamNames(key).head
    val reader: SamReader = SamReaderFactory.makeDefault().open(new File(filename))
    val it = reader.iterator
    var reads = new ArrayBuffer[Map[String, String]]
    while (it.hasNext) {
      val read = it.next
      reads += Map(
        "chrom" -> read.getReferenceName,
        "name" -> read.getReadName,
        "start" -> read.getAlignmentStart.toString,
        "end" -> read.getAlignmentEnd.toString,
        "cigar" -> read.getCigarString
        // ...
      )
    }
    Ok(toJson(reads))
  }


  /*
   * Return a URL pointing to a symbolic link to the requested portion of the BAM.
   * The function first extracts the region with "samtools view", writes the result to
   *   the local resource/temp, then creates a symbolic link to it in APACHE_TEMP_BAM_DIR.
   */
  def symlink(key: String, region: Option[String], description: Option[String]) = Action {

    def randomBamName(): String = "_" + (Random.alphanumeric take 19).mkString + ".bam"

    Logger.info(s"/symlink/$key?region=$region")
    val filename = getBamNames(key).head
    val bamOriginal: Path = Paths.get(BAM_PATH, filename)
    val randomName = randomBamName()

    cleanupOldTempFiles(TEMP_BAM_DIR,  regex = BAM_BAI_REGEX)
    cleanupOldTempFiles(APACHE_TEMP_BAM_DIR, regex = BAM_BAI_REGEX)

    if (! samtoolsExists()) InternalServerError("Could not find 'samtools' in $PATH.")
    else if (filename.isEmpty) InternalServerError(s"No corresponding BAM file for key '$key'")
    else if (! bamOriginal.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else {

    region match {
      // If no region, just add a symlink to the file where Apache can read it
      case None => Files.createSymbolicLink(
        Paths.get(APACHE_TEMP_BAM_DIR, randomName),
        Paths.get(BAM_PATH, filename)
      )
      // If a region is specified, extract the sub-BAM, write it to a temp directory, then create a symlink to it
      case Some(s: String) =>
        val dest: String = Paths.get(TEMP_BAM_DIR, randomName).toString
        val commandExtract = s"samtools view -hb $bamOriginal $s"  #>> new File(dest)
        Logger.info("Extracting region: " + commandExtract.toString)
        commandExtract.!
        indexBam(dest)
        Files.createSymbolicLink(
          Paths.get(APACHE_TEMP_BAM_DIR, randomName),
          Paths.get(TEMP_BAM_DIR, randomName)
        )
        Files.createSymbolicLink(
          Paths.get(APACHE_TEMP_BAM_DIR, randomName+".bai"),
          Paths.get(TEMP_BAM_DIR, randomName+".bai")
        )
    }

    Logger.info("Creating symlink to " + Paths.get(APACHE_TEMP_BAM_DIR, randomName).toString)
    val url = Paths.get(APACHE_TEMP_BAM_URL, randomName).toString
    Ok(url)
  }}


  /*
   * Download the file directly
   */
  //def download(key: String, region: Option[String], description: Option[String]) = Action {
  def download(key: String, region: String, description: Option[String]) = Action {  request =>
    Logger.info("---------------------------------------------------")
    Logger.info(s"/download/$key?region=$region")
    Logger.info("---------------------------------------------------")

    val bamRegex = """^([\w\d:-]+)(.bam){0,1}(.bai|.idx){0,1}$""".r
    val queryKey:String = key match {
      case bamRegex(root,b,null) => root
      case bamRegex(root,b,ext) => root
    }
    val queryRegion:String = region match {
      case "" => ""
      case bamRegex(root,b,null) => root
      case bamRegex(root,b,ext) => root
    }

    val filename = getBamNames(queryKey).head
    val bam: Path = Paths.get(BAM_PATH, filename)
    val index: Path = Paths.get(BAM_PATH, filename+".bai")
    val sha:String = hash(queryKey+queryRegion, "SHA")

    cleanupOldTempFiles(TEMP_BAM_DIR, regex = BAM_BAI_REGEX)
    if (! samtoolsExists()) InternalServerError("Could not find 'samtools' in $PATH.")
    else if (filename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bam.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else if (! index.toFile.exists) InternalServerError(s"BAM index not found on disk")
    else {

      region match {
        // If no region, send the whole file
        case "" =>
          Logger.info(">>>>>>>>  Returning full BAM")
          Ok.sendFile(bam.toFile)

        // If a region is specified, extract the sub-BAM, write it to a temp directory, and serve it
        //case Some(s: String) =>
        case s: String =>
          Logger.info(">>>>>>>>  Extract, index and send BAM/bai")
          val dest: Path = Paths.get(TEMP_BAM_DIR, sha+".bam")
          val idx: Path = Paths.get(TEMP_BAM_DIR, sha+".bam.bai")
          if (! dest.toFile.exists) {
            val commandExtract = s"samtools view -hb ${bam.toString} $queryRegion"  #>> dest.toFile
            Logger.info("Extracting region: " + commandExtract.toString)
            commandExtract.!
          }
          if (! idx.toFile.exists) {
            val commandIndex = s"samtools index ${dest.toString}"
            Logger.info("Indexing: " + commandIndex.toString)
            commandIndex.!
          }
          if (s.endsWith(".idx") || s.endsWith(".bai")) {
            Logger.info("BAI")
            Ok.sendFile(idx.toFile, fileName = _ => sha+"bam.bai")
          }
          else {
            Logger.info("BAM")
            Ok.sendFile(dest.toFile, fileName = _ => sha+"bam")
          }
      }

    }
  }


  def downloadIndex(key: String, description: Option[String]) = Action {
    Logger.info("---------------------------------------------------")
    Logger.info(s"/downloadIndex/$key")
    Logger.info("---------------------------------------------------")
    val regex = """^([\w\d:-]+)(.bam){0,1}(.bai|.idx){0,1}$""".r
    val queryKey: String = key match {
      case regex(root,b,ext) => root
      case _ => key
    }
    val bamFilename: String = getBamNames(queryKey).head
    val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
    val indexPath: Path = Paths.get(BAM_PATH, bamFilename+".bai")

    if (bamFilename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bamPath.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else {
      val indexExists = indexPath.toFile.exists
      if ((! indexExists) && (! samtoolsExists)) InternalServerError("Could not find 'samtools' in $PATH.")
      else {
        if (! indexExists) {
          indexBam(indexPath.toString)
        }
        Ok.sendFile(indexPath.toFile)
      }
    }
  }


  def downloadRange(key: String, description: Option[String]) = Action {  request =>
    Logger.info("---------------------------------------------------")
    Logger.info(s"/downloadRange/$key")
    Logger.info(request.toString)
    Logger.info("---------------------------------------------------")

    val bamRegex = """^([\w\d:-]+)(.bam){0,1}(.bai|.idx){0,1}$""".r
    val queryKey:String = key match {
      case bamRegex(root,b,ext) => root
      case _ => key
    }

    val bamFilename = getBamNames(queryKey).head
    val bam: File = Paths.get(BAM_PATH, bamFilename).toFile
    val index: File = Paths.get(BAM_PATH, bamFilename+".bai").toFile

    if (bamFilename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bam.exists) InternalServerError(s"BAM file not found on disk")
    else if (! index.exists) InternalServerError(s"BAM index not found on disk")
    else {

      val bamLength = bam.length
      val range = request.headers.get(RANGE) match {
        case None => None
        case Some(s:String) => s.replaceAll("bytes=", "").split("-").toList match {
          case rangeStart :: rangeEnd :: Nil => Some(rangeStart.toLong, rangeEnd.toLong)
          case rangeStart :: Nil => Some(rangeStart.toLong, bamLength)
          case _ => None
        }
      }
      val (start:Long, end:Long) = range.getOrElse(0L, bamLength)

      Logger.info(">>>>>>>>  Returning RANGE")
      assert(end-start < 500000, s"Query interval ($start-$end) too big")

      val status: Int = if (start != 0 || end != bamLength - 1) PARTIAL_CONTENT else OK
      val contentLength = if (status == PARTIAL_CONTENT) (end - start + 1) else bamLength

      val stream: FileInputStream = new FileInputStream(bam)
      stream.skip(start)
      val source = StreamConverters.fromInputStream(() => stream)

      val headers = mutable.Map(
        CONTENT_TYPE -> "application/octet-stream",
        ACCEPT_RANGES -> "bytes",
        CONTENT_LENGTH -> contentLength.toString,
        CONNECTION -> "keep-alive"
      )
      if (status == PARTIAL_CONTENT) {
        headers += CONTENT_RANGE -> s"bytes $start-$end/$contentLength"
      }
      Result(
        header = ResponseHeader(status, headers.toMap),
        body = HttpEntity.Streamed(source, Some(contentLength), Some("application/octet-stream"))
      )

    }
  }


}
