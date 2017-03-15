package controllers

import javax.inject.{Inject, _}

import controllers.generic.BamQueryController
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc._

import scala.util.{Failure, Success}


/**
 * Return reads from a BAM file in JSON format, i.e. an array of objects,
 * one object per aligned read.
 */
@Singleton
class JsonController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  //------------------ Actions -------------------//

  def bamPost(region: Option[String]) = Action { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        Ok("")
    }
  }

  def bamGet(sampleKey: String, token: String, region: Option[String]) = Action { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        Ok("")
    }
  }

}
