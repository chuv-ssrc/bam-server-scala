package controllers.bam

import javax.inject.{Inject, _}

import auth.AuthenticatedAction
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

  def AuthenticatedAction = new AuthenticatedAction(db)

  //------------------ Actions -------------------//

  def bamPost() = AuthenticatedAction(parse.json) { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        RangeResult.ofFile(br.bamFile, request.headers.get(RANGE), Some(BINARY))
    }
  }

  def bamGet(sample: String, token: Option[String], range: Option[String]) = AuthenticatedAction { implicit request =>
    sampleNameToBamRequest(sample) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        range match {
          // Maybe the range is in http headers
          case None => RangeResult.ofFile(br.bamFile, request.headers.get(RANGE), Some(BINARY))
          // If not, see if there is one in the url
          case Some(rangeFromUrl) =>
            val rangeRegex = """(\d*-\d*)""".r
            rangeRegex.findFirstIn(rangeFromUrl) match {
              case None => InternalServerError("Wrong format for argument 'range': should match /\\d*-\\d*/")
              case Some(r: String) =>
                RangeResult.ofFile(br.bamFile, Some("bytes="+r), Some(BINARY))
            }
        }
    }
  }

}
