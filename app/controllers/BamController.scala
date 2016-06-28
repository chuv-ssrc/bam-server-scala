package controllers

import javax.inject.Inject

import play.api.mvc._
import play.api.libs.json.Json._
import play.api.db._
import play.api.libs.iteratee.Enumerator
import htsjdk.samtools.{SAMFileHeader, SAMFileWriter, SAMFileWriterFactory, SamReader, SamReaderFactory}
import java.io.{ByteArrayInputStream, File}
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths}

import akka.stream.scaladsl.StreamConverters
import akka.util.ByteString
import play.api.http.HttpEntity
import play.api.libs.streams.Streams

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import sys.process._


class BamReader(val bam:String){
  val reader:SamReader = SamReaderFactory.makeDefault().open(new File(bam))

  def head() = {
    var i = 0
    val iterate = reader.iterator()
    while(iterate.hasNext && i < 10){
      System.err.println(iterate.next())
      i += 1
    }
  }
}


/**
  * Created by jdelafon on 06/03/16.
  */
class BamController @Inject()(db: Database) extends Controller {

  val bamPath = "resources/"
  val bamFilename = "141F.recal.bam"
  val testbam = bamPath + bamFilename

  /* Return the file names scorresponding to this key */
  def getBamNames(key: String) : List[String] = {
    db.withConnection { conn =>
      val cursor = conn.createStatement
      val res = cursor.executeQuery("SELECT `key`,`filename` FROM bam WHERE `key`='"+key+"';")
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      println(filenames)
      filenames.toList
    }
  }

  /* Return the bam text as output by "samtools view -hb <filename> <region>" */
  def read(key: String, region: String) = Action {
    val filename = getBamNames(key).head
    val bam = Paths.get(bamPath, filename)
    val command = s"samtools view -hb $bam $region" #| s"base64"
    val res = command.!!
    val bres = res.toCharArray.map(_.toByte)
    //println(res.length, res.isInstanceOf[Array[Byte]])
    // TODO: check the format of $region

    //Result(
    //  header = ResponseHeader(200, Map.empty),
    //  body = HttpEntity.Streamed(ByteString(res), Some(res.length), Some("text/plain"))
    //)

    Result(
      header = ResponseHeader(200, Map.empty),
      body = HttpEntity.Strict(ByteString(res), Some("application/octet-stream"))
    ).as(HTML)

    //Ok(res).as(HTML)

    //val reader: SamReader = SamReaderFactory.makeDefault().open(new File(testbam))
    //println(reader)
    ////val it = reader.iterator
    ////while (it.hasNext) { println(it.next()) }
    //val q = reader.query("1", 0, 7512448, true)
    //println(q)
    //reader.close()
    //val header: SAMFileHeader = reader.getFileHeader()
    //println(1, header)
    //val sf: SAMFileWriterFactory = new SAMFileWriterFactory
    //println(2, header, System.out, new File("/dev/stdout"))
    //val bw: SAMFileWriter = sf.makeBAMWriter(header, false, new File("/dev/stdout"))
    //println(3)
    //while (q.hasNext) {
    //  bw.addAlignment(q.next())
    //}

  }


  def view(key: String, region: String) = Action {
    val bam = new BamReader(testbam)
    val it = bam.reader.iterator
    var reads = new ArrayBuffer[Map[String, String]]
    while (it.hasNext) {
      val read = it.next
      reads += Map(
        "chrom" -> read.getReferenceName,
        "name" -> read.getReadName,
        "start" -> read.getAlignmentStart.toString,
        "end" -> read.getAlignmentEnd.toString,
        "cigar" -> read.getCigarString
      )
    }
    Ok(toJson(reads))
  }


  def symlink(key: String) = Action {
    val filename = getBamNames(key).head
    Files.createSymbolicLink(
      Paths.get("/Library/WebServer/Documents/bam", bamFilename),
      Paths.get(bamPath, bamFilename)
    )
    Ok("Creating symlink to " + Paths.get("/Library/WebServer/Documents/bam", bamFilename))
  }


  def download(key: String) = Action {
    val filename = getBamNames(key).head
    if (key == "asdf") {
      Ok.sendFile(
        content = new java.io.File(bamFilename),
        fileName = _ => bamFilename
      )
    } else {
      Ok("Got request" + key)
    }
  }

}
