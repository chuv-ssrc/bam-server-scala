package controllers

import java.io.File
import java.nio.file.{Path, Paths}
import javax.inject.Inject

import play.api.Logger
import play.api.db.Database
import play.api.mvc.{Action, RangeResult}
import utils.BamUtils._
import utils.Constants._
import utils.Utils._

import sys.process._



class BamDownloadController @Inject()(db: Database) extends BamController(db) {

  /**
    * Download the BAM index
    */
  def downloadIndex(key: String, description: Option[String]) = Action { request =>
    Logger.info("---------------------------------------------------")
    Logger.info("Request: " + request.toString)
    Logger.info("---------------------------------------------------")
    val queryKey: String = withoutBamExtension(key)
    val bamFilename: String = getBamName(queryKey)
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
        //Ok.sendFile(indexPath.toFile, inline=true)
        RangeResult.ofFile(indexPath.toFile, request.headers.get(RANGE), Some(BINARY))
      }
    }
  }


  /**
    * Respond to a query with a Range header.
    * Return the queried bytes range with a status 206 (Partial Content) and a ContentRange header.
    * If the query range is the whole file, return a status 200.
    *
    * An example of client sending such Range queries is IGV.js .
    * N.B. It is the job of the client to query or not more data if it receives 206.
    */
  def downloadRange(key: String, description: Option[String]) = Action {  request =>
    Logger.info("---------------------------------------------------")
    Logger.info("Request: " + request.toString)
    Logger.info("---------------------------------------------------")

    val queryKey: String = withoutBamExtension(key)
    val bamFilename: String = getBamName(queryKey)
    val bam: File = Paths.get(BAM_PATH, bamFilename).toFile

    if (bamFilename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database")
    else if (! bam.exists) InternalServerError(s"BAM file not found on disk")
    else {
      RangeResult.ofFile(bam, request.headers.get(RANGE), Some(BINARY))
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