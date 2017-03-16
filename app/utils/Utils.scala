package utils

import java.io.File
import java.security.MessageDigest

import play.api.Logger
import play.api.db.Database
import play.api.libs.json.JsValue
import play.api.mvc.{AnyContent, Request}

import scala.collection.mutable.ArrayBuffer
import sys.process._


object Utils {

  def stringValueFromRequestBody(request: Request[AnyContent], key: String): Option[String] = {
    request.body.asJson flatMap { body: JsValue =>
       (body \ key).asOpt[String]
    }
  }

  def hash(s: String, method: String = "SHA") = {
    MessageDigest.getInstance(method).digest((s).getBytes)
      .map(0xFF & _).map { "%02x".format(_) }.foldLeft(""){_ + _}
  }

  /**
    * Use the script onDisk.sh to check if the file exists,
    * because it is able to deal with long-term storage and mounted volumes like we have with /archive.
    */
  def isOnDisk(file: File, archive: Boolean = false): Boolean = {
    if (archive) {
      val isFound: Int = s"scripts/onDisk.sh ${file.toString}".!
      isFound == 0
    } else {
      file.exists
    }
  }

}


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
