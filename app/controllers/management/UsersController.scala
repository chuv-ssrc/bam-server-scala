package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedRequest}
import forms.UsersForm
import models.User
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.util.Try


@Singleton
class UsersController @Inject()(db: Database) extends Controller {

  val AdminAction = new AdminAction(db)

  /**
    * From the request body, make a list of Users - for insertion.
    * Give them the same appId as the admin user executing the action.
    *
    * N.B. In any case a web interface will be proper to a given app,
    * and a global admin can have an account by several auth providers.
    * There is no case when we want to add several users from multiple apps in one query.
    */
  def usersFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[User] = {
    // Get the app ID from the user executing the action;
    val appId = request.user.appId

    // Get the users part from the JSON body of the request
    val usersJs: JsArray = (request.body \ "users").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast users array from request body to JsArray")
    }
    // Transform the JSON body into a list of Users, with the same appId as the admin requesting this
    // - unless the appId is found in the JSON info for a given user.
    val users: Seq[User] = usersJs.value map { userJs =>
      Try (UsersForm.fromJson(userJs, appId)) getOrElse {
        throw new IllegalArgumentException("Could not cast user from request body to User model")
      }
    }
    users
  }

  /**
    * Extract the list of user names passed in request body, to a list of Strings - for deletion.
    */
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
    * Return the number of users found with this *username* and *appId* in the database.
    */
  def findUserByUsername(username: String, appId: Int): Int = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id` FROM `users` WHERE username = ? AND app_id = ? ;")
      statement.setString(1, username)
      statement.setInt(2, appId)
      val result = statement.executeQuery()
      var count = 0
      while (result.next()) {
        count += 1
      }
      count
    }
  }

  /**
    * Add one or more users to the database.
    * Expects a JSON body of the type
    * ```{ "users": [{"username": "A"}, {"username": "B"}] }```
    */
  def addUsers() = AdminAction(parse.json) { implicit request =>

    val users = usersFromRequest
    // they all get an appId from the admin requesting this

    db.withConnection { conn =>
      val unknowns: String = users.map(_ => "(?,?)").mkString(",")
      val statement = conn.prepareStatement("INSERT INTO `users`(`app_id`,`username`) VALUES "+unknowns+" ;")
      users.zipWithIndex.foreach { case (user, i: Int) =>
        statement.setInt(2*i+1, user.appId)  // user.appId is always defined in this case
        statement.setString(2*i+2, user.username)
      }
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

    db.withConnection { conn =>
      val unknowns: String = usernames.map(_ => "?").mkString(",")
      val statement = conn.prepareStatement("DELETE FROM `users` WHERE `username` IN ("+unknowns+") ;")
      usernames.zipWithIndex.foreach {case (username, i: Int) => statement.setString(i+1, username)}
      statement.execute()
      Ok(s"Deleted ${usernames.size} user(s)")
    }
  }

}
