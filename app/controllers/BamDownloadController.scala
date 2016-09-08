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

    if (bamFilename.isEmpty) InternalServerError(s"No corresponding BAM file for that key in database.")
    else {
      val indexFile = indexPath.toFile
      if (! isOnDisk(indexFile)) {
        if (! samtoolsExists())
          InternalServerError("Index file not found, and could not find 'samtools' in $PATH to fix it.")
        else if (! isOnDisk(bamPath.toFile)) {
          InternalServerError("Index file not found, and BAM file not found on disk.")
        }
        else {
          indexBam(indexPath.toString)
          RangeResult.ofFile(indexFile, request.headers.get(RANGE), Some(BINARY))
        }
      } else {
        RangeResult.ofFile(indexFile, request.headers.get(RANGE), Some(BINARY))
      }
    }
  }

  /**
    * Download the BAM file. Respond to a query with a Range header.
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
    else if (! isOnDisk(bam)) InternalServerError(s"BAM file not found on disk")
    else {
      RangeResult.ofFile(bam, request.headers.get(RANGE), Some(BINARY))
    }
  }


}