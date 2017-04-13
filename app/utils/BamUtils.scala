package utils

import play.api.Logger
import play.api.db.Database
import scala.collection.mutable.ArrayBuffer
import sys.process._


/**
  * Operations related to BAM files, such as inv
  */
object BamUtils {

  val stdout = new StringBuilder
  val stderr = new StringBuilder
  val logger = ProcessLogger(stdout append _, stderr append _)

  /**
    * Return whether samtools is found in PATH
    * @return Boolean
    */
  def samtoolsExists(): Boolean = {
    val exists = "which samtools" ! logger  == 0
    if (! exists) {
      Logger.error("Could not find 'samtools' in $PATH.")
      Logger.error(System.getenv("PATH"))
    }
    exists
  }

  /**
    * Run samtools index on the give bam file.
    * @param bamFilename: BAM file name
    */
  def indexBam(bamFilename:String): Unit = {
    val commandIndex = s"samtools index $bamFilename"
    Logger.info("Indexing: " + commandIndex.toString)
    commandIndex ! logger
  }

  /**
    * Return the file name corresponding to this key in the database
    */
  def getBamName(db: Database, key: String): String = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `name`,`filename` FROM `samples` WHERE `name`=?;")
      statement.setString(1, key)
      val res = statement.executeQuery()
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      //Logger.debug(s"File names for key $key: ${filenames.mkString}")
      if (filenames.isEmpty) "" else filenames.head
    }
  }


}
