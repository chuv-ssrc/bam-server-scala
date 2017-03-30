package controllers.generic

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc.{AnyContent, Controller, Request}
import utils.BamUtils.getBamName
import utils.Common.isOnDisk
import utils.Common.stringValueFromRequestBody
import scala.util.{Failure, Success, Try}



class BamQueryController @Inject()(db: Database, config: Configuration) extends Controller {

  val BAM_PATH = config.getString("env.BAM_PATH").get
  require(Paths.get(BAM_PATH).toFile.exists, s"BAM_PATH not found ($BAM_PATH). Check your configuration file (conf/application.conf).")

  /**
    * Try to get the JSON body of a POST request;
    * Try to get the "key" argument from the body;
    * Make a BamRequest from that - see `keyToBamRequest`.
    */
  def parseBamRequestFromPost(request: Request[AnyContent]): Try[BamRequest] = {
    stringValueFromRequestBody(request, "key") match {
      case None =>
        Failure(new Exception(s"No key found in request body."))
      case Some(key: String) =>
        keyToBamRequest(key)
    }
  }

  /**
    * Check that the sample key is found in the database;
    * Check that the file it refers to exists/is readable on disk;
    * If so, return a Success(BamRequest).
    */
  def keyToBamRequest(sampleKey: String): Try[BamRequest] = {
    val bamFilename: String = getBamName(db, sampleKey)
    val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
    val indexPath: Path = Paths.get(BAM_PATH, bamFilename + ".bai")
    if (bamFilename.isEmpty) {
      Failure(new Exception(s"No corresponding BAM file for key $sampleKey in database."))
    } else if (!isOnDisk(bamPath.toFile)) {
    //} else if (!bamPath.toFile.exists) {
      Failure(new Exception(s"This BAM file cannot be found at BAM_PATH=$BAM_PATH."))
    } else {
      Success(BamRequest(bamPath.toFile, indexPath.toFile))
    }
  }


}