package controllers

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import htsjdk.samtools.{SamReader, SamReaderFactory}
import play.api.db.Database
import play.api.libs.json.Json._
import play.api.mvc.Action
import utils.Constants._

import scala.collection.mutable.ArrayBuffer


/**
  * Return a JSON summarizing the content of this BAM region (1 item per alignment)
  */
class BamViewController @Inject()(db: Database) extends BamController(db) {

  def view(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamName(key)
    val bamPath: Path = Paths.get(BAM_PATH, filename)
    val reader: SamReader = SamReaderFactory.makeDefault().open(bamPath.toFile)
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

}