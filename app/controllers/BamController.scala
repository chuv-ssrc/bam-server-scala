package controllers

import java.io.File
import java.nio.file.Paths
import javax.inject.Inject

import play.api.mvc._
import play.api.db._
import play.api.Logger
import utils.Constants._

import sys.process._

import scala.collection.mutable.ArrayBuffer


/**
  * Generic controller for BAM queries
  * @param db: database name
  */
class BamController @Inject()(db: Database) extends Controller {

  /**
    * Return the file name corresponding to this key in the database
    */
  def getBamName(key: String): String = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `key`,`filename` FROM bam WHERE `key`=?;")
      statement.setString(1, key)
      val res = statement.executeQuery()
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      Logger.debug(s"File names for key $key: ${filenames.mkString}")
      if (filenames.isEmpty) "" else filenames.head
    }
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

