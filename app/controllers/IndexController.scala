package controllers

import java.nio.file.{Path, Paths}
import javax.inject._
import java.io.File
import java.nio.file.{Path, Paths}
import javax.inject.Inject

import play.api.{Configuration, Logger}
import play.api.db.Database
import play.api.mvc.{Action, Controller, RangeResult}
import play.api.libs.json._
import utils.BamUtils._
import utils.Utils._
import scala.Predef.require


/**
 * This controller creates an `Action` to handle HTTP requests to the
 * application's home page.
 */
@Singleton
class IndexController @Inject()(db: Database, config: Configuration) extends Controller {

  val BAM_PATH = config.getString("env.BAM_PATH").get
  require(Paths.get(BAM_PATH).toFile.exists, "BAM_PATH not found. Exiting.")


  def baiPost = Action { implicit request =>
    //Logger.debug(request.toString)
    val token = request.headers.get("Authorization")
    request.body.asJson match {
      case None =>
        InternalServerError("No body was found in POST request")
      case Some(body: JsValue) =>
        (body \ "key").asOpt[String] match {
          case None =>
            InternalServerError(s"No key found in request body.")
          case Some(key: String) =>
            val bamFilename: String = getBamName(db, key)
            val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
            val indexPath: Path = Paths.get(BAM_PATH, bamFilename + ".bai")
            if (bamFilename.isEmpty) {
              //Logger.error(s"Key '$key' not found in database")
              InternalServerError(s"No corresponding BAM file for that key in database.")
            } else {
              val indexFile = indexPath.toFile
              if (!isOnDisk(indexFile)) {
                if (!samtoolsExists())
                  InternalServerError("Index file not found, and could not find 'samtools' in $PATH to fix it.")
                else if (!isOnDisk(bamPath.toFile)) {
                  InternalServerError("Index file not found, and BAM file not found on disk.")
                }
                else {
                  indexBam(indexPath.toString)
                  RangeResult.ofFile(indexFile, request.headers.get(RANGE), Some(BINARY))
                }
              } else {
                RangeResult.ofFile(indexFile, request.headers.get(RANGE), Some(BINARY))
              }
            }
        }
    }

  }

  def baiGet(token: String) = Action {
    //Ok(views.html.index("BAM server operational."))
    Ok("BAM server operational.")
  }

}
