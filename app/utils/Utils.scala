package utils

import java.io.File
import java.security.MessageDigest
import java.util.Calendar

import play.api.Logger
import play.api.db.Database

import scala.collection.mutable.ArrayBuffer
import sys.process._


object Utils {

  def hash(s:String, method:String="SHA") = {
    MessageDigest.getInstance(method).digest((s).getBytes)
      .map(0xFF & _).map { "%02x".format(_) }.foldLeft(""){_ + _}
  }

  /** Delete all files in *dir* matching this *regex* that are older than *seconds*
   * Default is one hour, all extensions.
   */
  def cleanupOldTempFiles(dir:String, seconds:Long=3600, regex:String=".*"): Unit = {
    val now:Long = Calendar.getInstance().getTime.getTime
    val delta:Long = seconds * 1000
    //Logger.info(s"Erasing all files in '$dir' older than ${delta / 1000 / 60} min")
    new File(dir).listFiles.filter(_.getName.matches(".*?\\.bam.*"))
      .filter(_.lastModified < now - delta).foreach(_.delete)
  }

  /**
   * Remove '.bam' or '.bam.bai' extensions from the input string
   */
  def withoutBamExtension(s: String) : String = {
    val queryKey:String = s match {
      case Constants.KEY_REGEX(root,b,ext) => root
      case _ => s
    }
    queryKey
  }

  def isOnDisk(file: File, archive: Boolean = false): Boolean = {
    if (archive) {
      val isFound: Int = s"scripts/onDisk.sh ${file.toString}".!
      isFound == 0
    } else (
      file.exists
      )
  }

}


object BamUtils {

  /**
    * Return whether samtools is found in PATH
    * @return Boolean
    */
  def samtoolsExists(): Boolean = {
    val exists = "which samtools".! == 0
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
    commandIndex.!
  }

  /**
    * Return the file name corresponding to this key in the database
    */
  def getBamName(db: Database, key: String): String = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `key`,`filename` FROM bam WHERE `key`=?;")
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
