package controllers.bam

import javax.inject._

import auth.AuthenticatedAction
import controllers.generic.BamQueryController
import scala.concurrent.ExecutionContext.Implicits.global
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc._
import utils.BamUtils._
import utils.Common._

import scala.concurrent.Future
import scala.util.{Failure, Success}


/**
 * Provide the BAM index file (.bai).
 */
@Singleton
class IndexController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  def AuthenticatedAction = new AuthenticatedAction(db)

  //------------------ Actions -------------------//

  def baiPost = AuthenticatedAction.async { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        // Logger.debug(err.getMessage)
        Future(InternalServerError(err.getMessage))
      case Success(br: BamRequest) =>
        getBamIndex(br)
    }
  }

  def baiGet(sample: String, token: Option[String]) = AuthenticatedAction.async { implicit request =>
    keyToBamRequest(sample) match {
      case Failure(err) =>
        // Logger.debug(err.getMessage)
        Future(InternalServerError(err.getMessage))
      case Success(br: BamRequest) =>
        getBamIndex(br)
    }
  }

  //--------------- Common code ----------------//

  /**
    * Async because of `indexBam`.
    */
  def getBamIndex(br: BamRequest)(implicit request: Request[AnyContent]): Future[Result] = {
    /* Index not found, try to index */
    if (!isOnDisk(br.indexFile)) {
      /* Cannot index because no samtools not found */
      if (!samtoolsExists()) {
        Future(InternalServerError("Index file not found, and could not find 'samtools' in $PATH to fix it."))
        /* Index the bam using samtools */
      } else {
        indexBam(br.indexFile.toPath.toString) recover {
          case err => InternalServerError("Error while indexing: "+ err.getMessage)
        } map { _ =>
          Ok.sendFile(br.indexFile)
        }
      }
    /* Index found, return */
    } else {
      Future.successful(Ok.sendFile(br.indexFile))
    }
  }

}
