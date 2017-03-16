package controllers

import javax.inject.{Inject, _}

import controllers.generic.BamQueryController
import htsjdk.samtools.{SamReader, SamReaderFactory}
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.collection.mutable.ArrayBuffer
import scala.util.{Failure, Success}


/**
  * Return reads from a BAM file in JSON format, i.e. an array of objects,
  * one object per aligned read.
  * The returned fields correspond to the BAM file columns:
  * - name: read name (QNAME)
  * - flag: (FLAG)
  * - chrom: (RNAME)
  * - start: leftmost aligned position; 1-based, 0 if unmapped (POS)
  * - end: rightmost aligned position, i.e. start + read length; 1-based, 0 if unmapped
  * - mapq: mapping quality (MAPQ)
  * - cigar: cigar string (CIGAR)
  * - rnext: ? (RNEXT)
  * - pnext: ? (PNEXT)
  * - tlen: "template length", aka "insert size" (TLEN)
  * - seq: read sequence (SEQ)
  * - qual: quality string (per-base) (QUAL)
  */
@Singleton
class JsonController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  implicit val anyWrites = new Writes[Any] {
    def writes(x: Any) = x match {
      case n: Int => JsNumber(n)
      case s: String => JsString(s)
      case b: Boolean => JsBoolean(b)
      case m: Map[String, Any] @unchecked => JsObject(m.mapValues(writes))
      case a: Seq[Any] => JsArray(a.map(writes))
      case _ => JsNull
    }
  }

  //------------------ Actions -------------------//

  def bamPost(region: Option[String]) = Action { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        readsToJson(br, region)
    }
  }

  def bamGet(sampleKey: String, token: String, region: Option[String]) = Action { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        readsToJson(br, region)
    }
  }

  //--------------- Common code ----------------//

  /**
    * See https://samtools.github.io/htsjdk/javadoc/htsjdk
    * under SamRecord for all available attributes.
    */
  def readsToJson(br: BamRequest, region: Option[String]) = {
    val reader: SamReader = SamReaderFactory.makeDefault().open(br.bamFile)
    val it = reader.iterator
    var reads = new ArrayBuffer[Map[String, Any]]
    while (it.hasNext) {
      val read = it.next
      reads += Map(
        "name" -> read.getReadName,
        "flag" -> read.getFlags,
        "chrom" -> read.getReferenceName,
        "start" -> read.getAlignmentStart,
        "end" -> read.getAlignmentEnd,  // start + read length
        "mapq" -> read.getMappingQuality,
        "cigar" -> read.getCigarString,
        //"rnext" -> read.getMateReferenceName,
        //"pnext" -> read.getMateReferenceIndex.toInt,
        "tlen" -> read.getInferredInsertSize,  // "insert size", aka "template length"
        "seq" -> read.getReadString,
        "qual" -> read.getBaseQualityString
      )
    }
    Ok(Json.toJson(reads))
  }

}
