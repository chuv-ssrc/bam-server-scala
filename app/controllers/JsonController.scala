package controllers

import javax.inject.{Inject, _}

import auth.AuthenticatedAction
import controllers.generic.BamQueryController
import htsjdk.samtools.{SAMRecord, SamReader, SamReaderFactory}
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.collection.mutable.ArrayBuffer
import scala.util.{Failure, Success}
import utils.JsonUtils.anyWrites


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

  //------------------ Actions -------------------//

  def bamPost(region: Option[String]) = AuthenticatedAction { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        readsToJson(br, region)
    }
  }

  def bamGet(sampleKey: String, token: Option[String], region: Option[String]) = AuthenticatedAction { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        readsToJson(br, region)
    }
  }

  //--------------- Common code ----------------//

  def mapFromRead(read: SAMRecord) = {
    Map(
      "name" -> read.getReadName,
      "flag" -> read.getFlags,
      "chrom" -> read.getReferenceName,
      "start" -> read.getAlignmentStart,
      "end" -> read.getAlignmentEnd, // start + read length
      "mapq" -> read.getMappingQuality,
      "cigar" -> read.getCigarString,
      "rnext" -> (if (read.getMateReferenceName == read.getReferenceName) "=" else read.getMateReferenceName),
      "pnext" -> read.getMateAlignmentStart,
      "tlen" -> read.getInferredInsertSize, // "insert size", aka "template length"
      "seq" -> read.getReadString,
      "qual" -> read.getBaseQualityString
    )
  }

  def parseRegion(region: String): (String, Int, Int) = {
    val regex = """(\w+):([0-9]+)-([0-9]+)""".r
    region match {
      case regex(chrom: String, start: String, end: String) => (chrom, start.toInt, end.toInt)
      case _ => throw new IllegalArgumentException(s"Invalid region: '$region'")
    }
  }

  def isInRange(read: SAMRecord, region: Option[String]): Boolean = {
    region match {
      case None => true
      case Some(reg: String) =>
        val (chrom, start, end) = parseRegion(reg)
        read.getReferenceName == chrom && read.getAlignmentStart >= start && read.getAlignmentEnd <= end
      case _ => false
    }
  }

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
      if (isInRange(read, region)) {
        reads += mapFromRead(read)
      }
    }
    Ok(Json.toJson(reads))
  }

}
