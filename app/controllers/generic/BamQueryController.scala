package controllers.generic

import java.nio.file.{Path, Paths}
import javax.inject.Inject

import auth.AuthenticatedRequest
import models.{BamRequest, User}
import play.api.Configuration
import play.api.db.Database
import play.api.libs.json.JsValue
import play.api.mvc.Controller
import utils.BamUtils.getBamName
import utils.Common.isOnDisk

import scala.collection.mutable.ListBuffer
import scala.util.{Failure, Success, Try}



class BamQueryController @Inject()(db: Database, config: Configuration) extends Controller {

  val BAM_PATH = config.getString("env.BAM_PATH").get
  require(Paths.get(BAM_PATH).toFile.exists, s"BAM_PATH not found ($BAM_PATH). Check your configuration file (conf/application.conf).")

  /**
    * Try to get the JSON body of a POST request;
    * Try to get the "sample" argument from the body;
    * Make a BamRequest from that - see `sampleNameToBamRequest` below.
    */
  def parseBamRequestFromPost(request: AuthenticatedRequest[JsValue]): Try[BamRequest] = {
    (request.body \ "sample").asOpt[String] match {
      case None =>
        Failure(new Exception(s"No key found in request body."))
      case Some(sampleName: String) =>
        val user: User = request.user
        val canAccess: Boolean = canUserAccessSample(user, sampleName)
        sampleNameToBamRequest(sampleName)
    }
  }

  /**
    * Check that the sample key is found in the database;
    * Check that the file it refers to exists/is readable on disk;
    * If so, return a Success(BamRequest).
    */
  def sampleNameToBamRequest(sampleName: String): Try[BamRequest] = {
    val bamFilename: String = getBamName(db, sampleName)
    val bamPath: Path = Paths.get(BAM_PATH, bamFilename)
    val indexPath: Path = Paths.get(BAM_PATH, bamFilename + ".bai")
    if (bamFilename.isEmpty) {
      Failure(new Exception(s"No corresponding BAM file for key $sampleName in database."))
    } else if (!isOnDisk(bamPath.toFile)) {
    //} else if (!bamPath.toFile.exists) {
      Failure(new Exception(s"This BAM file cannot be found at BAM_PATH=$BAM_PATH."))
    } else {
      Success(BamRequest(bamPath.toFile, indexPath.toFile))
    }
  }

  /**
    * Return whether this *user* has the right to see this sample.
    * @param user: User making a request for some sample.
    * @param sampleName: name of the sample the user tries to access.
    */
  def canUserAccessSample(user: User, sampleName: String): Boolean = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement(s"""
        SELECT count(*) FROM `users_samples`
          JOIN `samples` ON samples.id = users_samples.sample_id
          JOIN `users` ON users.id = users_samples.user_id
        WHERE samples.name = ?
          AND users.username = ?
          AND users.app_id = ?
          AND users_samples.isActive = 1 AND users.isActive = 1 AND samples.isActive = 1
      """)
      statement.setString(1, sampleName)
      statement.setString(2, user.username)
      statement.setInt(3, user.appId)
      val res = statement.executeQuery()
      var accesses = 0
      while (res.next()) {
        accesses += res.getInt(1)
      }
      accesses > 0
    }
  }

}



