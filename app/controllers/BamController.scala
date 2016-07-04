package controllers

import java.io.File
import java.util.Calendar
import java.security.MessageDigest
import javax.inject.Inject

import play.api.mvc._
import play.api.libs.json.Json._
import play.api.db._
import play.api.Logger
import play.api.Environment._
import play.api.libs.iteratee.Enumerator
import play.api.http.HttpEntity
//import play.api.libs.streams.Streams
//import akka.util.ByteString

import java.io.{ByteArrayInputStream, File}
import java.nio.file.{Files, Path, Paths}
//import java.nio.charset.StandardCharsets

//import scala.io.Source
//import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.util.Random

import htsjdk.samtools.{SamReader, SamReaderFactory}
//import htsjdk.samtools.{SAMFileHeader, SAMFileWriter, SAMFileWriterFactory}
import com.typesafe.config.ConfigFactory
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

  val BAM_PATH:String = ConfigFactory.load().getString("env.BAM_PATH")
  val APACHE_TEMP_BAM_DIR:String = ConfigFactory.load().getString("env.APACHE_TEMP_BAM_DIR")
  val APACHE_TEMP_BAM_URL:String = ConfigFactory.load().getString("env.APACHE_TEMP_BAM_URL")
  val TEMP_BAM_DIR:String = Paths.get(ConfigFactory.load().getString("env.TEMP_BAM_DIR")).toAbsolutePath.toString

  def samtoolsExists(): Boolean = {
    val exists = "which samtools".! == 0
    if (! exists) {
      Logger.error("Could not find 'samtools' in $PATH.")
      Logger.error(System.getenv("PATH"))
    }
    exists
  }

  def hash(s:String, method:String="SHA") = {
    MessageDigest.getInstance(method).digest((s).getBytes)
      .map(0xFF & _).map { "%02x".format(_) }.foldLeft(""){_ + _}
  }

  def randomBamName(): String = "_" + (Random.alphanumeric take 19).mkString + ".bam"

  /* Return the file names scorresponding to this key */
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

  def cleanupTempBams(dir:String): Unit = {
    val now:Long = Calendar.getInstance().getTime.getTime
    val HOUR:Long = 3600 * 1000
    val delta:Long = 1 * HOUR  // number of milliseconds before deleting
    //Logger.info(s"Erasing all files in '$dir' older than ${delta / 1000 / 60} min")
    new File(dir).listFiles.filter(_.getName.matches(".*?\\.bam.*"))
      .filter(_.lastModified < now - delta).foreach(_.delete)
  }

  /*
   * Return the BAM text as output by "samtools view -hb <filename> <region>".
   * It is encoded in base64 because transfers have to be ascii (?).
   */
  def read(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamNames(key).head
    val bam: Path = Paths.get(BAM_PATH, filename)
    val command = s"samtools view -hb $bam $region" #| "base64"
    val res: String = command.!!
    // TODO: check the format of $region
    Ok(res).as(HTML)

    //val stream = new ByteArrayInputStream(res.getBytes)
    //Ok.stream(Enumerator.fromStream(stream))
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

  def indexBam(bamFilename:String): Unit = {
    val commandIndex = s"samtools index $bamFilename"
    Logger.info("Indexing: " + commandIndex.toString)
    commandIndex.!
  }

  /*
   * Return a URL pointing to a symbolic link to the requested portion of the BAM.
   * The function first extracts the region with "samtools view", writes the result to
   *   the local resource/temp, then creates a symbolic link to it in APACHE_TEMP_BAM_DIR.
   */
  def symlink(key: String, region: Option[String], description: Option[String]) = Action {
    Logger.info(s"/symlink/$key?region=$region")
    val filename = getBamNames(key).head
    val bamOriginal: Path = Paths.get(BAM_PATH, filename)
    val randomName = randomBamName()

    cleanupTempBams(TEMP_BAM_DIR)
    cleanupTempBams(APACHE_TEMP_BAM_DIR)

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


  def sendFile(file:File, name:String="dummy") = {
    assert(file.exists, s"File not found: ${file.getAbsolutePath.toString}")
    Ok.sendFile(
      content = file,
      fileName = _ => name,
      inline = true
    )
      //.as(HTML)
      .withHeaders(
      CONNECTION -> "keep-alive",
      ACCEPT_RANGES -> "bytes",
      ACCEPT_ENCODING -> "gzip, deflate, sdch",
      ACCEPT_LANGUAGE -> "en-US,en;q=0.8,fr;q=0.6,ru;q=0.4,de;q=0.2,pt;q=0.2"
    )
  }

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
      //case None => ""
      //case Some(bamRegex(root,null)) => root
      //case Some(bamRegex(root,ext)) => root
      case "" => ""
      case bamRegex(root,b,null) => root
      case bamRegex(root,b,ext) => root
    }

    val filename = getBamNames(queryKey).head
    val bam: Path = Paths.get(BAM_PATH, filename)
    val index: Path = Paths.get(BAM_PATH, filename+".bai")
    val sha:String = hash(queryKey+queryRegion, "SHA")
    println(BAM_PATH)
    println(queryKey, queryRegion, sha)

    cleanupTempBams(TEMP_BAM_DIR)
    if (! samtoolsExists()) InternalServerError("Could not find 'samtools' in $PATH.")
    else if (filename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bam.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else if (! index.toFile.exists) InternalServerError(s"BAM index not found on disk")
    else {


      val range = request.headers.get(RANGE) match {
        case None => None
        case Some(s:String) => s.replaceAll("bytes=", "").split("-").toList match {
          case rangeStart :: rangeEnd :: Nil => Some(rangeStart.toLong, rangeEnd.toLong)
          case rangeStart :: Nil => Some(rangeStart.toLong, bam.toFile.length)
          case _ => None
        }
      }


      region match {
        // If no region, send the whole file
        //case None =>
        case "" =>
          Logger.info(">>>>>>>>  Returning full BAM")
          sendFile(bam.toFile, sha+".bam")

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
          if (s.endsWith(".idx") | s.endsWith(".bai")) {
            Logger.info("BAI")
            sendFile(idx.toFile, sha+".bam.bai")
          }
          else {
            Logger.info("BAM")
            sendFile(dest.toFile, sha+".bam")
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
      case regex(root,b, null) => root
      case _ => key
    }
    val bamFilename: String = getBamNames(queryKey).head
    val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
    val indexPath: Path = Paths.get(BAM_PATH, bamFilename+".bai")

    if (bamFilename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bamPath.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else {
      val indexExists = indexPath.toFile.exists
      if ((! indexExists) & (! samtoolsExists)) InternalServerError("Could not find 'samtools' in $PATH.")
      else {
        if (! indexExists) {
          indexBam(indexPath.toString)
        }

        sendFile(indexPath.toFile)

      }




    }

    Ok("asdf")
  }


  def downloadRange(key: String, description: Option[String]) = Action {  request =>
    Logger.info("---------------------------------------------------")
    Logger.info(s"/downloadRange/$key")
    Logger.info("---------------------------------------------------")

    val bamRegex = """^([\w\d:-]+)(.bam){0,1}(.bai|.idx){0,1}$""".r
    val queryKey:String = key match {
      case bamRegex(root,b,null) => root
      case bamRegex(root,b,ext) => root
    }

    val filename = getBamNames(queryKey).head
    val bam: Path = Paths.get(BAM_PATH, filename)
    val index: Path = Paths.get(BAM_PATH, filename+".bai")

    cleanupTempBams(TEMP_BAM_DIR)
    if (! samtoolsExists()) InternalServerError("Could not find 'samtools' in $PATH.")
    else if (filename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bam.toFile.exists) InternalServerError(s"BAM file not found on disk")
    else if (! index.toFile.exists) InternalServerError(s"BAM index not found on disk")
    else {

      val range = request.headers.get(RANGE) match {
        case None => None
        case Some(s:String) => s.replaceAll("bytes=", "").split("-").toList match {
          case rangeStart :: rangeEnd :: Nil => Some(rangeStart.toLong, rangeEnd.toLong)
          case rangeStart :: Nil => Some(rangeStart.toLong, bam.toFile.length)
          case _ => None
        }
      }

      Logger.info(">>>>>>>>  Returning RANGE")
      sendFile(bam, filename)

    }
  }


}
