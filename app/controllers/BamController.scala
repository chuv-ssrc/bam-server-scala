package controllers

import javax.inject.Inject

import play.api.mvc._
import play.api.db._
import play.api.Logger

import scala.collection.mutable.ArrayBuffer


/**
  * Generic controller for BAM queries
  * @param db: database name
  */
class BamController @Inject()(db: Database) extends Controller {

  /**
    * Return the file names scorresponding to this key in the database
    */
  def getBamName(key: String): String = {
    db.withConnection { conn =>
      val cursor = conn.createStatement
      val res = cursor.executeQuery("SELECT `key`,`filename` FROM bam WHERE `key`='" + key + "';")
      val filenames = ArrayBuffer.empty[String]
      while (res.next()) {
        filenames += res.getString("filename")
      }
      Logger.debug(s"File names for key $key: ${filenames.mkString}")
      if (filenames.isEmpty) "" else filenames.head
    }
  }

}




