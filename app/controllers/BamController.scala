package controllers

import javax.inject.Inject

import play.api.mvc._
import play.api.libs.json.Json._
import play.api.db._
import play.api.http.HttpEntity
import play.api.libs.streams.Streams
import akka.util.ByteString

import java.io.{ByteArrayInputStream, File}
import java.nio.file.{Files, Path, Paths}
import java.nio.charset.StandardCharsets

import scala.io.Source
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.util.Random

import htsjdk.samtools.{SAMFileHeader, SAMFileWriter, SAMFileWriterFactory, SamReader, SamReaderFactory}
import com.typesafe.config.ConfigFactory
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

  val BAM_PATH = ConfigFactory.load().getString("env.BAM_PATH")
  val APACHE_BAM_DIR = ConfigFactory.load().getString("env.APACHE_BAM_DIR")
  val TEMP_BAM_DIR = ConfigFactory.load().getString("env.TEMP_BAM_DIR")

  /* Return the file names scorresponding to this key */
  def getBamNames(key: String) : List[String] = {
    db.withConnection { conn =>
      val cursor = conn.createStatement
      val res = cursor.executeQuery("SELECT `key`,`filename` FROM bam WHERE `key`='"+key+"';")
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      println("File names: " + filenames.mkString)
      filenames.toList
    }
  }

  /*
   * Return the BAM text as output by "samtools view -hb <filename> <region>".
   * It is encoded in base64 because transfers have to be ascii (?).
   */
  def read(key: String, region: String) = Action {
    val filename = getBamNames(key).head
    val bam = Paths.get(BAM_PATH, filename)
    val command = s"samtools view -hb $bam $region" #| "base64"
    val res = command.!!
    // TODO: check the format of $region
    Ok(res).as(HTML)
  }


  /*
   * Return a JSON summarizing the content of this BAM region (1 item per alignment)
   */
  def view(key: String, region: String) = Action {
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
   *   the local resource/temp, then creates a symbolic link to it in APACHE_BAM_DIR.
   */
  def symlink(key: String, region: Option[String]) = Action {
    val filename = getBamNames(key).head
    val bamOriginal = Paths.get(BAM_PATH, filename)
    val randomName = "_" + (Random.alphanumeric take 19).mkString + ".bam"

    region match {
      // If no region, just add a symlink to the file where Apache can read it
      case None => Files.createSymbolicLink(
        Paths.get(APACHE_BAM_DIR, randomName),
        Paths.get(BAM_PATH, filename)
      )
      // If a region is specified, extract the sub-BAM, write it to a temp directory, then create a symlink to it
      case Some(s: String) =>
        val dest: String = Paths.get(TEMP_BAM_DIR, randomName).toString
        val commandExtract = s"samtools view -hb $bamOriginal $s"  #>> new File(dest)
        val commandIndex = s"samtools index $dest"
        println("Extracting region: " + commandExtract.toString)
        println("Indexing: " + commandIndex.toString)
        commandExtract.!
        commandIndex.!
        Files.createSymbolicLink(
          Paths.get(APACHE_BAM_DIR, randomName),
          Paths.get(TEMP_BAM_DIR, randomName)
        )
        Files.createSymbolicLink(
          Paths.get(APACHE_BAM_DIR, randomName+".bai"),
          Paths.get(TEMP_BAM_DIR, randomName+".bai")
        )
    }

    println("Creating symlink to " + Paths.get(APACHE_BAM_DIR, randomName).toString)
    val url = "http://localhost/bam/" + randomName
    Ok(url)
  }


  /*
   * Download the whole file directly
   */
  def download(key: String, region: String) = Action {
    val filename = getBamNames(key).head
    val bam = Paths.get(BAM_PATH, filename).toString
    Ok.sendFile(
      content = new java.io.File(bam),
      fileName = _ => filename
    )
  }

}
