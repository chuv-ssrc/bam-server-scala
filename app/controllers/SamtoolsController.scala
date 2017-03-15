package controllers

import javax.inject.{Inject, _}
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
 */
@Singleton
class SamtoolsController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  //------------------ Actions -------------------//

  def bamPost(maybeRegion: Option[String]) = Action { implicit request =>
    parseBamRequestFromPost(request) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
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

  def bamGet(sampleKey: String, token: String, region: Option[String]) = Action { implicit request =>
    keyToBamRequest(sampleKey) match {
      case Failure(err) =>
        //Logger.debug(err.getMessage)
        InternalServerError(err.getMessage)
      case Success(br: BamRequest) =>
        Ok("")
    }
  }


//    /**
//      * Return the BAM text as output by "samtools view -hb <filename> <region>".
//      * It is encoded in base64 because transfers have to be ascii (?).
//      * Subject to "java.lang.OutOfMemoryError: Java heap space" when the region is too big.
//      */
//    def read(key: String, region: String, description: Option[String]) = Action {
//      val filename = getBamName(key)
//      val bamPath: Path = Paths.get(BAM_PATH, filename)
//      if (! samtoolsExists())
//        InternalServerError("Could not find 'samtools' in $PATH.")
//      else {
//        val command = s"samtools view -hb $bamPath $region" #| "base64"
//        val res: String = command.!!
//        Ok(res).as(HTML)
//      }
//    }

}
