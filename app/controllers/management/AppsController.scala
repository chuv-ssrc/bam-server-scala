package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedRequest}
import forms.AppForm
import models.App
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.util.Try


@Singleton
class AppsController @Inject()(db: Database) extends Controller {

  val AdminAction = new AdminAction(db)

  def appFromRequest(implicit request: AuthenticatedRequest[JsValue]): App = {
    val app: App =
      Try (AppForm.fromJson(request.body)) getOrElse {
        throw new IllegalArgumentException("Could not cast request body to App model")
      }
    app
  }

  def issFromRequest(implicit request: AuthenticatedRequest[JsValue]): String = {
    val iss: String = (request.body \ "iss").asOpt[String] getOrElse {
      throw new IllegalArgumentException("Could not find iss claim in request body")
    }
    iss
  }

  def findAppByIss(iss: String): Int = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id` FROM `apps` WHERE iss = ? ;")
      statement.setString(1, iss)
      val result = statement.executeQuery()
      var count = 0
      while (result.next()) {
        count += 1
      }
      count
    }
  }

  /**
    * Add an app to the database.
    * Expects a JSON body of the type
    * ```{ "iss": "as in jwt claim", "keyFile": "/", "description": "" }```
    */
  def addApp() = AdminAction(parse.json) { implicit request =>

    val app: App = appFromRequest

    if (findAppByIss(app.iss) > 0) {
      InternalServerError(s"Cannot insert app '${app.iss}' because it already exists")
    } else {
      db.withConnection { conn =>
        val statement = conn.prepareStatement("INSERT INTO `apps`(`iss`,`keyFile`,`description`) VALUES (?,?,?) ;")
        statement.setString(1, app.iss)
        statement.setString(2, app.keyFile)
        statement.setString(3, app.description.getOrElse(""))
        statement.execute()
        Ok(s"Inserted app '${app.iss}'")
      }
    }
  }

  /**
    * Delete an app from the database.
    * Expects a JSON body of the type
    * ```{ "iss": "as in jwt claim" }```
    */
  def deleteApp() = AdminAction(parse.json) { implicit request =>

    val iss = issFromRequest

    if (findAppByIss(iss) == 0) {
      InternalServerError(s"Cannot delete app '$iss' because it does not exist")
    } else {
      db.withConnection { conn =>
        val statement = conn.prepareStatement("DELETE FROM `apps` WHERE `iss` = ? ;")
        statement.setString(1, iss)
        statement.execute()
        Ok(s"Deleted app '$iss'")
      }
    }
  }

}
