package controllers

import javax.inject.Inject

import play.api.mvc.{Action, Controller}
import play.api.libs.json.Json._
import play.api.db._
import htsjdk.samtools.{SamReader, SamReaderFactory}
import java.io.File

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer


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

  val bamFilename = "resources/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"

  def view = Action {
    val bam = new BamReader(bamFilename)
    val it = bam.reader.iterator
    var reads = new ArrayBuffer[Map[String, String]]
    while (it.hasNext) {
      val read = it.next
      reads += Map(
        "name" -> read.getReadName,
        "start" -> read.getAlignmentStart.toString,
        "end" -> read.getAlignmentEnd.toString,
        "cigar" -> read.getCigarString
      )
    }
    Ok(toJson(reads))
  }

  def download(name: String) = Action { implicit request =>
    db.withConnection { conn =>
      val cursor = conn.createStatement
      val res = cursor.executeQuery("SELECT name,filename from variants_db")
      val m = new mutable.HashMap[String, String]
      while (res.next()) {
        m += (res.getString("name") -> res.getString("filename"))
      }
      println(m)
    }

    if (name == "asdf") {
      Ok.sendFile(
        content = new java.io.File(bamFilename),
        fileName = _ => "HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam"
      )
    } else {
      Ok("Got request [" + request + "] " + name)
    }
  }

}
