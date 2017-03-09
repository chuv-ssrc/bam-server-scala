package controllers

import java.nio.file.{Path, Paths}
import javax.inject.{Inject, _}

import controllers.generic.BamQueryController
import play.api.Configuration
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc.{Action, Controller, RangeResult}
import utils.BamUtils._
import utils.Utils._

import scala.Predef.require


/**
 * This controller creates an `Action` to handle HTTP requests to the
 * application's home page.
 */
@Singleton
class RangeController @Inject()(db: Database, config: Configuration) extends BamQueryController(db, config) {

  def bamPost(range: String) = Action { implicit request =>
      Ok{""}
  }

  def bamGet(sampleKey: String, range: String, token: String) = Action {
    Ok("BAM server operational.")
  }

}
