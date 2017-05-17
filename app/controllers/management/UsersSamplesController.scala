package controllers.management

import javax.inject._

import auth.{AdminAction, AuthenticatedAction, AuthenticatedRequest}
import forms.UsersSamplesForm
import models.UserSample
import play.api.db.Database
import play.api.libs.json._
import play.api.mvc._
import scala.util.Try


@Singleton
class UsersSamplesController @Inject()(db: Database) extends Controller {

  val AdminAction = new AdminAction(db)
  import UsersSamplesController._

  /**
    * Add one or more users_samples rows to the database.
    * Expects a JSON body of the type
    * ```{ "users": [{"username": "A"}, {"username": "B"}] }```
    */
  def addUsersSamples() = AdminAction(parse.json) { implicit request =>

    val appId = request.user.appId
    val usersSamples = usersSamplesFromRequest
    val userIds = getUserIds(db, usersSamples.map(_.username))
    val sampleIds = getSampleIds(db, usersSamples.map(_.sample))

    val counts: Seq[Int] = (userIds zip sampleIds).map(us => findUsersSamples(db, us._1, us._2))
    val notExistYet = (userIds zip sampleIds).zipWithIndex.filter(x => counts(x._2) == 0).map(_._1).zipWithIndex

    if (notExistYet.nonEmpty) {
      db.withConnection { conn =>
        val unknowns: String = notExistYet.map(_ => "(?,?)").mkString(",")
        val statement = conn.prepareStatement(s"""
          INSERT INTO `users_samples`(`user_id`,`sample_id`) VALUES $unknowns ;
        """)
        notExistYet.foreach { case ((uid, sid), i) =>
          statement.setInt(2 * i + 1, uid)
          statement.setInt(2 * i + 2, sid)
        }
        statement.execute()
      }
    }
    Ok(s"Inserted ${notExistYet.size} user-sample(s)")
  }

  /**
    * Delete one or more users_samples rows from the database.
    * Expects a JSON body of the type
    * ```{ "users": ["A","B"] }```
    */
  def deleteUsersSamples() = AdminAction(parse.json) { implicit request =>

    val appId = request.user.appId
    val usersSamples = usersSamplesFromRequest
    val userIds = getUserIds(db, usersSamples.map(_.username))
    val sampleIds = getSampleIds(db, usersSamples.map(_.sample))

    val counts: Seq[Int] = (userIds zip sampleIds).map(us => findUsersSamples(db, us._1, us._2))
    val alreadyExist = (userIds zip sampleIds).zipWithIndex.filter(x => counts(x._2) > 0)
    // It is no problem in SQL to try and delete things that don't exist

    db.withConnection { conn =>
      (userIds zip sampleIds).foreach { case (uid, sid) =>
        val statement = conn.prepareStatement(s"""
          DELETE FROM `users_samples` WHERE user_id = ? and sample_id = ? ;
        """)
        statement.setInt(1, uid)
        statement.setInt(2, sid)
        statement.execute()
      }
      Ok(s"Deleted ${alreadyExist.size} user-sample(s)")
    }
  }

}



object UsersSamplesController {

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

  def getUserIds(db: Database, usernames: Seq[String]): Seq[Int] = {
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

  def getSampleIds(db: Database, sampleNames: Seq[String]): Seq[Int] = {
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
    * Return the number of users found with this *username* and *appId* in the database.
    */
  def findUsersSamples(db: Database, userId: Int, sampleId: Int): Int = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id` FROM `users_samples` WHERE user_id = ? AND sample_id = ? ;")
      statement.setInt(1, userId)
      statement.setInt(2, sampleId)
      val result = statement.executeQuery()
      var count = 0
      while (result.next()) {
        count += 1
      }
      count
    }
  }

}

