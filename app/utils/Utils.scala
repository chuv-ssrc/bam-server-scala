package utils

import java.io.File
import java.security.MessageDigest
import java.util.Calendar

import play.api.Logger

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

}
