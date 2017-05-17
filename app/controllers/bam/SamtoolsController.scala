package controllers.bam

import javax.inject.{Inject, _}

import auth.AuthenticatedAction
import controllers.generic.BamQueryController
import scala.concurrent.ExecutionContext.Implicits.global
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc.Result
import utils.BamUtils.samtoolsExists

import scala.concurrent.Future
import scala.sys.process._
import scala.util.{Failure, Success}


/**
 * Provide a slice of a BAM file as returned by `samtools view -hb <bam> <region>`
 * in the body of the response.
 * If the given region is out of range, return an empty bam.
 */
@Singleton
class SamtoolsController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  def AuthenticatedAction = new AuthenticatedAction(db)

  //------------------ Actions -------------------//

  def bamPost(region: Option[String]) = AuthenticatedAction.async(parse.json) { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err: IllegalAccessException) =>
        Future(Unauthorized(err.getMessage))
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        Future(InternalServerError(err.getMessage))
      case Success(br: BamRequest) =>
        getBamWithSamtools(br, region)
    }
  }

  def bamGet(sample: String, token: Option[String], region: Option[String]) = AuthenticatedAction.async { implicit request =>
    sampleNameToBamRequest(sample, request.user) match {
      case Failure(err: IllegalAccessException) =>
        Future(Unauthorized(err.getMessage))
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        Future(InternalServerError(err.getMessage))
      case Success(br: BamRequest) =>
        getBamWithSamtools(br, region)
    }
  }

  //--------------- Common code ----------------//

  def getBamWithSamtools(br: BamRequest, maybeRegion: Option[String]): Future[Result] = {
    maybeRegion match {
      case None =>
        Future.successful(Ok.sendFile(br.bamFile))
      case Some(region) =>
        if (! samtoolsExists())
          Future(InternalServerError("Could not find 'samtools' in $PATH."))
        else {
          Future {
            val command = s"samtools view -hb ${br.bamFile.toPath} $region"
            val res: String = command.!!
            Ok(res).as(BINARY)
          }
        }
    }
  }

}
