package controllers

import javax.inject.{Inject, _}

import auth.AuthenticatedAction
import controllers.generic.BamQueryController
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc._
import utils.BamUtils.samtoolsExists

import scala.util.{Failure, Success}
import sys.process._


/**
 * Provide a slice of a BAM file as returned by `samtools view -hb <bam> <region>`
 * in the body of the response.
 * If the given region is out of range, return an empty bam.
 */
@Singleton
class SamtoolsController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  //------------------ Actions -------------------//

  def bamPost(region: Option[String]) = AuthenticatedAction { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamWithSamtools(br, region)
    }
  }

  def bamGet(sampleKey: String, token: Option[String], region: Option[String]) = AuthenticatedAction { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamWithSamtools(br, region)
    }
  }

  //--------------- Common code ----------------//

  def getBamWithSamtools(br: BamRequest, maybeRegion: Option[String]) = {
    maybeRegion match {
      case None =>
        Ok.sendFile(br.bamFile)
      case Some(region) =>
        if (! samtoolsExists())
          InternalServerError("Could not find 'samtools' in $PATH.")
        else {
          val command = s"samtools view -hb ${br.bamFile.toPath} $region"
          val res: String = command.!!
          Ok(res).as(BINARY)
        }
    }
  }

}
