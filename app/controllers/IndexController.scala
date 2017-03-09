package controllers

import javax.inject._

import controllers.generic.BamQueryController
import models.BamRequest
import play.api.{Configuration, Logger}
import play.api.db.Database
import play.api.mvc._
import utils.BamUtils._
import utils.Utils._

import scala.util.{Failure, Success}


/**
 * This controller creates an `Action` to handle HTTP requests to the
 * application's home page.
 */
@Singleton
class IndexController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  def baiPost = Action { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamIndex(br)
    }

  }

  def baiGet(sampleKey: String, token: String) = Action { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamIndex(br)
    }
  }

  //--------------- Common code ----------------//

  def getBamIndex(br: BamRequest)(implicit request: Request[AnyContent]): Result = {
    //Logger.debug(request.toString)
    val token = request.headers.get("Authorization")
    /* Index not found, try to index */
    if (!isOnDisk(br.indexFile)) {
      /* Cannot index because no samtools not found */
      if (!samtoolsExists()) {
        InternalServerError("Index file not found, and could not find 'samtools' in $PATH to fix it.")
        /* Index the bam using samtools */
      } else {
        indexBam(br.indexFile.toPath.toString)
        RangeResult.ofFile(br.indexFile, request.headers.get(RANGE), Some(BINARY))
      }
      /* Index found, return */
    } else {
      RangeResult.ofFile(br.indexFile, request.headers.get(RANGE), Some(BINARY))
    }
  }

}
