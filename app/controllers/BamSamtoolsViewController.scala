package controllers

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import play.api.Logger
import play.api.db.Database
import play.api.mvc.Action
import utils.Constants._
import utils.BamUtils.samtoolsExists
import utils.Utils._

import sys.process._


class BamSamtoolsViewController @Inject()(db: Database) extends BamController(db) {

  /**
    * Return the BAM text as output by "samtools view -hb <filename> <region>".
    * It is encoded in base64 because transfers have to be ascii (?).
    * Subject to "java.lang.OutOfMemoryError: Java heap space" when the region is too big.
    */
  def read(key: String, region: String, description: Option[String]) = Action {
    val filename = getBamName(key)
    val bamPath: Path = Paths.get(BAM_PATH, filename)
    if (! samtoolsExists())
      InternalServerError("Could not find 'samtools' in $PATH.")
    else {
      val command = s"samtools view -hb $bamPath $region" #| "base64"
      val res: String = command.!!
      Ok(res).as(HTML)
    }
  }

  /**
    * Download the part of the BAM file corresponding to this region.
    * Writes temporarily the extracted region to `env.TEMP_BAM_DIR`.
    * Old "classic" version that does it with samtools.
    */
  def download(key: String, region: String, description: Option[String]) = Action {  request =>
    Logger.info("---------------------------------------------------")
    Logger.info(s"/download/$key?region=$region")
    Logger.info("---------------------------------------------------")

    val queryKey:String = key match {
      case KEY_REGEX(root,b,null) => root
      case KEY_REGEX(root,b,ext) => root
    }
    val queryRegion:String = region match {
      case "" => ""
      case KEY_REGEX(root,b,null) => root
      case KEY_REGEX(root,b,ext) => root
    }

    val filename = getBamName(queryKey)
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

}

