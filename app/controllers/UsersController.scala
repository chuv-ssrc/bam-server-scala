package controllers

import javax.inject._

import models.User
import play.api.db.Database
import play.api.mvc._
import auth.AuthenticatedAction
import play.api.libs.json._
import forms.UsersForm

import scala.util.Try


@Singleton
class UsersController @Inject()(db: Database) extends Controller {

  val AuthenticatedAction = new AuthenticatedAction(db)

  /**
    * Add one or more users to the database.
    */
  def addUsers() = AuthenticatedAction(parse.json) { request =>

    // Get the app ID from the user executing the action;
    //  he cannot add anything on behalf of another app anyway.
    val appId = request.user.appId

    val usersJs: JsArray = (request.body \ "users").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast users array from request body to JsArray")
    }
    val users: Seq[User] = usersJs.value map { userJs =>
      Try (UsersForm.fromJson(userJs, appId)) getOrElse {
        throw new IllegalArgumentException("Could not cast user from request body to User model")
      }
    }

    db.withTransaction { conn =>
      for { user <- users } yield {
        val statement = conn.prepareStatement("INSERT INTO `users`(`username`) VALUES (?);")
        statement.setString(1, user.username)
        statement.execute()
      }
      Ok(s"Inserted ${users.size} users")
    }
  }

}
