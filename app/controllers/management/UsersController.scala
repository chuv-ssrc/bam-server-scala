package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedAction, AuthenticatedRequest}
import forms.UsersForm
import models.User
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.util.Try


@Singleton
class UsersController @Inject()(db: Database) extends Controller {

  val AuthenticatedAction = new AuthenticatedAction(db)
  val AdminAction = new AdminAction(db)

  def usersFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[User] = {
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
    users
  }

  def usernamesFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[String] = {
    val usernamesJs: JsArray = (request.body \ "users").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast usernames array from request body to JsArray")
    }
    val usernames: Seq[String] = usernamesJs.value map { usernameJs =>
      Try (usernameJs.as[String]) getOrElse {
        throw new IllegalArgumentException("Could not cast username from request body to String")
      }
    }
    usernames
  }

  /**
    * Add one or more users to the database.
    * Expects a JSON body of the type
    * ```{ "users": [{"username": "A"}, {"username": "B"}] }```
    */
  def addUsers() = AdminAction(parse.json) { implicit request =>

    val users = usersFromRequest

    db.withTransaction { conn =>
      val unknowns: String = users.map(_ => "(?)").mkString(",")
      val statement = conn.prepareStatement("INSERT INTO `users`(`username`) VALUES "+unknowns+" ;")
      users.zipWithIndex.foreach {case (user, i: Int) => statement.setString(i+1, user.username)}
      statement.execute()
      Ok(s"Inserted ${users.size} user(s)")
    }
  }

  /**
    * Delete one or more users from the database.
    * Expects a JSON body of the type
    * ```{ "users": ["A","B"] }```
    */
  def deleteUsers() = AdminAction(parse.json) { implicit request =>

    val usernames: Seq[String] = usernamesFromRequest

    db.withTransaction { conn =>
      val unknowns: String = usernames.map(_ => "?").mkString(",")
      val statement = conn.prepareStatement("DELETE FROM `users` WHERE `username` IN ("+unknowns+") ;")
      usernames.zipWithIndex.foreach {case (username, i: Int) => statement.setString(i+1, username)}
      statement.execute()
      Ok(s"Deleted ${usernames.size} user(s)")
    }
  }

}
