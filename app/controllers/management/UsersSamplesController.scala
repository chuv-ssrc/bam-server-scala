package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedAction, AuthenticatedRequest}
import forms.UsersSamplesForm
import models.UserSample
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._

import scala.collection.mutable.ArrayBuffer
import scala.util.Try


@Singleton
class UsersSamplesController @Inject()(db: Database) extends Controller {

  val AuthenticatedAction = new AuthenticatedAction(db)
  val AdminAction = new AdminAction(db)

  def usersSamplesFromRequest(implicit request: AuthenticatedRequest[JsValue]): Seq[UserSample] = {
    val usersSamplesJs: JsArray = (request.body \ "users_samples").asOpt[JsArray] getOrElse {
      throw new IllegalArgumentException("Could not cast users_samples array from request body to JsArray")
    }
    val users: Seq[UserSample] = usersSamplesJs.value map { userSampleJs =>
      Try (UsersSamplesForm.fromJson(userSampleJs)) getOrElse {
        throw new IllegalArgumentException("Could not cast request body to UserSample models")
      }
    }
    users
  }

  def getUserIds(usernames: Seq[String]): Seq[Int] = {
    val unknowns: String = usernames.map(_ => "?").mkString(",")
    db.withConnection { conn =>
      val statement = conn.prepareStatement(s"SELECT username,id FROM users WHERE username IN ($unknowns) ;")
      usernames.zipWithIndex.foreach {case (name, i) => statement.setString(i+1, name)}
      val res = statement.executeQuery()
      val usersMap = scala.collection.mutable.HashMap[String,Int]()
      while (res.next()) {
        usersMap += (res.getString("username") -> res.getInt("id"))
      }
      val userIds = usernames map (usersMap(_))
      userIds
    }
  }

  def getSampleIds(sampleNames: Seq[String]): Seq[Int] = {
    val unknowns: String = sampleNames.map(_ => "?").mkString(",")
    db.withConnection { conn =>
      val statement = conn.prepareStatement(s"SELECT name,id FROM samples WHERE name IN ($unknowns) ;")
      sampleNames.zipWithIndex.foreach {case (name, i) => statement.setString(i+1, name)}
      val res = statement.executeQuery()
      val samplesMap = scala.collection.mutable.HashMap[String,Int]()
      while (res.next()) {
        samplesMap += (res.getString("name") -> res.getInt("id"))
      }
      val sampleIds = sampleNames map (samplesMap(_))
      sampleIds
    }
  }

  /**
    * Add one or more users_samples rows to the database.
    * Expects a JSON body of the type
    * ```{ "users": [{"username": "A"}, {"username": "B"}] }```
    */
  def addUsersSamples() = AdminAction(parse.json) { implicit request =>

    val usersSamples = usersSamplesFromRequest
    val userIds = getUserIds(usersSamples.map(_.username))
    val sampleIds = getSampleIds(usersSamples.map(_.sample))

    db.withTransaction { conn =>
      val unknowns: String = usersSamples.map(_ => "(?,?)").mkString(",")
      val statement = conn.prepareStatement(s"""
        INSERT INTO `users_samples`(`user_id`,`sample_id`) VALUES $unknowns ;
      """)
      (userIds zip sampleIds).zipWithIndex.foreach { case ((uid, sid), i) =>
        statement.setInt(2*i+1, uid)
        statement.setInt(2*i+2, sid)
      }
      statement.execute()
      Ok(s"Inserted ${usersSamples.size} users")
    }
  }

  /**
    * Delete one or more users_samples rows from the database.
    * Expects a JSON body of the type
    * ```{ "users": ["A","B"] }```
    */
  def deleteUsersSamples() = AdminAction(parse.json) { implicit request =>

    val usersSamples = usersSamplesFromRequest
    val userIds = getUserIds(usersSamples.map(_.username))
    val sampleIds = getSampleIds(usersSamples.map(_.sample))

    db.withTransaction { conn =>
      (userIds zip sampleIds).foreach { case (uid, sid) =>
        val statement = conn.prepareStatement(s"""
          DELETE FROM `users_samples` WHERE user_id = ? and sample_id = ? ;
        """)
        statement.setInt(1, uid)
        statement.setInt(2, sid)
        statement.execute()
      }
      Ok(s"Deleted ${usersSamples.size} users")
    }
  }

}
