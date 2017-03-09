package controllers.generic

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.libs.json.JsValue
import play.api.mvc.{AnyContent, Controller, Request}
import utils.BamUtils.getBamName
import utils.Utils.isOnDisk

import scala.util.{Failure, Success, Try}



class BamQueryController @Inject()(db: Database, config: Configuration) extends Controller {

  val BAM_PATH = config.getString("env.BAM_PATH").get
  require(Paths.get(BAM_PATH).toFile.exists, s"BAM_PATH not found ($BAM_PATH). Check your configuration file (conf/application.conf).")

  def parseKeyFromPostRequest(request: Request[AnyContent]): Try[String] = {
    request.body.asJson match {
      case None =>
        Failure(new Exception("No body was found in POST request"))
      case Some(body: JsValue) =>
        (body \ "key").asOpt[String] match {
          case None =>
            Failure(new Exception(s"No key found in request body."))
          case Some(key: String) =>
            Success(key)
        }
    }
  }

  def parseBamRequest(request: Request[AnyContent]): Try[BamRequest] = {
    parseKeyFromPostRequest(request) match {
      case Failure(err) => Failure(err)
      case Success(key: String) =>
        val bamFilename: String = getBamName(db, key)
        val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
        val indexPath: Path = Paths.get(BAM_PATH, bamFilename + ".bai")
        if (bamFilename.isEmpty) {
          Failure(new Exception(s"No corresponding BAM file for key $key in database."))
        } else if (!isOnDisk(bamPath.toFile)) {
          Failure(new Exception(s"This BAM file cannot be found at BAM_PATH=$BAM_PATH."))
        } else {
          Success(BamRequest(bamPath.toFile, indexPath.toFile))
        }
    }
  }

}