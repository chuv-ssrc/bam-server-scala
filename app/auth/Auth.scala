package auth

import java.io.File
import java.security.PublicKey
import javax.inject._

import pdi.jwt.{JwtAlgorithm, JwtJson, JwtOptions}

import scala.util.Try
import play.api.libs.json._
import models.User
import utils.Common._
import auth.RSA._
import play.api.db.Database


@Singleton
class Auth @Inject()(db: Database) {

  private def getPublicKey(keyFile: String, path: String = "resources/rsa_keys"): PublicKey = {
    /* Try to find a certificate in *path* with extension ".cer" */
    findInTree(path, "cer") find {_.getName == keyFile} map {
      cert: File => readPublicKeyFromCertificate(cert.getPath)
    } getOrElse {
      /* Try to find a public key in *path* with extension ".pem" */
      findInTree(path, "pem") find {_.getName == keyFile} map {
        keyFile: File => readPublicKeyFromPemFile(keyFile.getPath)
      } getOrElse {
        throw new IllegalArgumentException("Could not find a suitable RSA public key in either .cer or .pem file.")
      }
    }
  }

  private def issFromClaim(claim: JsObject): String = {
    (claim \ "iss").asOpt[String] getOrElse {
      throw new IllegalArgumentException("No 'iss' claim found in token (to identify the issuer)")
    }
  }

  private def userFromClaim(claim: JsObject): String = {
    (claim \ "name").asOpt[String] getOrElse {
      throw new IllegalArgumentException("No 'name' claim found in token (username)")
    }
  }

  private def appFromDatabase(iss: String): (Int, String) = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id`,`keyFile` FROM `apps` WHERE `iss`=?;")
      statement.setString(1, iss)
      val result = statement.executeQuery()
      var res: Option[(Int, String)] = None
      while (result.next()) {
        val appId = result.getInt("id")
        val keyFile = result.getString("keyFile")
        res = Some((appId, keyFile))
      }
      res getOrElse {
        throw new IllegalArgumentException("Could not find an app with this 'iss' identifer in database")
      }
    }
  }

  private def userFromDatabase(username: String, appId: Int): User = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT * FROM `users` WHERE `username`=? AND `app_id`=?;")
      statement.setString(1, username)
      statement.setInt(2, appId)
      val result = statement.executeQuery()
      var res: Option[User] = None
      while (result.next()) {
        val name = result.getString("username")
        val appid = result.getInt("app_id")
        val isActive = result.getBoolean("isActive")
        res = Some(User(name, Some(appid), isActive))
      }
      res getOrElse {
        throw new IllegalArgumentException("Could not find a user with this name and app_id in the database")
      }
    }
  }

  def validateToken(jwt: String, db: Database): Try[User] = {
    // Get the claim without verification, to read the app name from the "iss" claim
    // and read the corresponding public key, then validate the token
    JwtJson.decodeJson(jwt, JwtOptions(signature=false)) flatMap { claim =>
      for {
        iss <- Try(issFromClaim(claim))
        username <- Try(userFromClaim(claim))
        (appId, keyFile) <- Try(appFromDatabase(iss))
        user <- Try(userFromDatabase(username, appId))
        publicKey <- Try(getPublicKey(keyFile))
      } yield {
        // Now validate the token with the public key
        // Raises an error if the validation fails
        JwtJson.validate (jwt, publicKey, Seq(JwtAlgorithm.RS256))
        user
      }
    }
  }

}
