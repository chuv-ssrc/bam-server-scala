package controllers

import javax.inject.{Inject, _}
import controllers.generic.BamQueryController
import models.BamRequest
import play.api.Configuration
import play.api.db.Database
import play.api.mvc._
import scala.util.{Failure, Success}


/**
 * Provide a slice of a BAM file by reading a bytes range header, such as
 * "Range: bytes=0-65639", from the HTTP request.
 * If no Range header is present, the whole file is returned (200).
 * If the bytes range overlaps the file, a partial content is returned (206).
 * If the range does not overlap, an error is thrown (416).
 */
@Singleton
class RangeController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  //------------------ Actions -------------------//

  def bamPost() = Action { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamRange(br)
    }
  }

  def bamGet(sampleKey: String, range: String, token: String) = Action { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        getBamRange(br)
    }
  }

  //--------------- Common code ----------------//

  def getBamRange(br: BamRequest)(implicit request: Request[AnyContent]): Result = {
    RangeResult.ofFile(br.bamFile, request.headers.get(RANGE), Some(BINARY))
  }

}
