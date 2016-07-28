package controllers

import java.io.File
import java.nio.file.{Files, Path, Paths}
import javax.inject.Inject

import play.api.Logger
import play.api.db.Database
import play.api.mvc.Action
import utils.BamUtils._
import utils.Constants._
import utils.Utils._

import scala.util.Random

import sys.process._


/**
  * Return a URL pointing to a symbolic link to the requested portion of the BAM.
  * The function first extracts the region with "samtools view", writes the result to
  * `TEMP_BAM_DIR`, then creates a symbolic link to it in APACHE_TEMP_BAM_DIR`.
  */
class BamSymlinkController @Inject()(db: Database) extends BamController(db) {

  def symlink(key: String, region: Option[String], description: Option[String]) = Action {

    def randomBamName(): String = "_" + (Random.alphanumeric take 19).mkString + ".bam"

    Logger.info(s"/symlink/$key?region=$region")
    val filename = getBamName(key)
    val bamOriginal: Path = Paths.get(BAM_PATH, filename)
    val randomName = randomBamName()

    cleanupOldTempFiles(TEMP_BAM_DIR, regex = BAM_BAI_REGEX)
    cleanupOldTempFiles(APACHE_TEMP_BAM_DIR, regex = BAM_BAI_REGEX)

    if (!samtoolsExists()) InternalServerError("Could not find 'samtools' in $PATH.")
    else if (filename.isEmpty) InternalServerError(s"No corresponding BAM file for key '$key'")
    else if (!bamOriginal.toFile.exists) InternalServerError(s"BAM file not found on disk")
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
          val commandExtract = s"samtools view -hb $bamOriginal $s" #>> new File(dest)
          Logger.info("Extracting region: " + commandExtract.toString)
          commandExtract.!
          indexBam(dest)
          Files.createSymbolicLink(
            Paths.get(APACHE_TEMP_BAM_DIR, randomName),
            Paths.get(TEMP_BAM_DIR, randomName)
          )
          Files.createSymbolicLink(
            Paths.get(APACHE_TEMP_BAM_DIR, randomName + ".bai"),
            Paths.get(TEMP_BAM_DIR, randomName + ".bai")
          )
      }

      Logger.info("Creating symlink to " + Paths.get(APACHE_TEMP_BAM_DIR, randomName).toString)
      val url = Paths.get(APACHE_TEMP_BAM_URL, randomName).toString
      Ok(url)
    }
  }
}