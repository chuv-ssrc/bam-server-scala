package auth

import java.io.File
import java.security.PublicKey
import javax.inject._

import pdi.jwt.{JwtAlgorithm, JwtJson, JwtOptions}

import scala.util.Try
import play.api.libs.json._
import models.User
import utils.Common._
import utils.Constants._
import auth.RSA._
import play.api.db.Database


@Singleton
class Auth @Inject()(db: Database) {

  /*
  /**
    * Given the name of an RSA public key file, return the PublicKey object.
    * @param keyPath: the name of the file containing the public RSA key, including the .crt/.pem extension.
    * @param path: the path to the directory containing the RSA public key certificates.
    */
  private def getPublicKey(keyPath: String, path: String = "resources/rsa_keys"): PublicKey = {
    /* Try to find a certificate in *path* with extension ".cer" */
    findInTree(path, "cer") find {_.getName == keyPath} map { cert: File =>
      readPublicKeyFromCertificate(cert.getPath)
    } getOrElse {
      /* Try to find a public key in *path* with extension ".pem" */
      findInTree(path, "pem") find {_.getName == keyPath} map { keyFile: File =>
        readPublicKeyFromPemFile(keyFile.getPath)
      } getOrElse {
        throw new IllegalArgumentException("Could not find a suitable RSA public key in either .cer or .pem file.")
      }
    }
  }

  private def getPublicKeyFromDb(appId: Int): PublicKey = {
    readPublicKeyFromDb(db, appId)
  }
  */

  /**
    * Read the "iss" (issuer) claim from JWT and return its value (the app name).
    * @param claim
    */
  private def issFromClaim(claim: JsObject): String = {
    (claim \ TOKEN_ISS_KEY).asOpt[String] getOrElse {
      throw new IllegalArgumentException("No 'iss' claim found in token (to identify the issuer)")
    }
  }

  /**
    * Read the "username" claim from JWT and return its value (the user name).
    * @param claim
    */
  private def userFromClaim(claim: JsObject): String = {
    (claim \ TOKEN_USER_KEY).asOpt[String] getOrElse {
      throw new IllegalArgumentException("No 'name' claim found in token (username)")
    }
  }

  /**
    * Search the database for an app with identifier *iss*.
    * Return a couple (appId, key).
    * @param iss: the app name/identifier, as it is in the JWT under the "iss" claim.
    */
  private def appFromDatabase(iss: String): (Int, String, String) = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT `id`,`key`,`algorithm` FROM `apps` WHERE `iss`=?;")
      statement.setString(1, iss)
      val result = statement.executeQuery()
      var res: Option[(Int, String, String)] = None
      while (result.next()) {
        val appId = result.getInt("id")
        val key = result.getString("key")
        val algorithm = result.getString("algorithm")
        res = Some((appId, key, algorithm))
      }
      res getOrElse {
        throw new IllegalArgumentException(s"Could not find app '$iss' in database")
      }
    }
  }

  /**
    * Search the database for a user with identifier *username* for this *appId*.
    * Return a User object.
    */
  private def userFromDatabase(username: String, appId: Int): User = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("SELECT * FROM `users` WHERE `username`=? AND `app_id`=?;")
      statement.setString(1, username)
      statement.setInt(2, appId)
      val result = statement.executeQuery()
      var res: Option[User] = None
      while (result.next()) {
        val isActive = result.getBoolean("isActive")
        val isAdmin = result.getBoolean("isAdmin")
        res = Some(User(appId, username, isActive=isActive, isAdmin=isAdmin))
      }
      res getOrElse {
        throw new IllegalArgumentException(s"Could not find user '$username' with app_id=$appId in the database")
      }
    }
  }

  def validateToken(jwt: String, keyString: String, algorithm: String): Try[Unit] = {
    if (algorithm.startsWith("RS")) {
      Try(readPublicKeyFromString(keyString.replace("\\n","\n"))).map { key =>
        algorithm match {
          case "RS256" => JwtJson.validate(jwt, key, Seq(JwtAlgorithm.RS256))
          case "RS384" => JwtJson.validate(jwt, key, Seq(JwtAlgorithm.RS256))
          case "RS512" => JwtJson.validate(jwt, key, Seq(JwtAlgorithm.RS256))
          case _ => throw new IllegalArgumentException("Only algorithms RS256,RS384,RS512,HS256,HS384,HS512 are supported.")
        }
      }
    } else {
      Try(
        algorithm match {
          case "HS256" => JwtJson.validate(jwt, keyString, Seq(JwtAlgorithm.HS256))
          case "HS384" => JwtJson.validate(jwt, keyString, Seq(JwtAlgorithm.HS384))
          case "HS512" => JwtJson.validate(jwt, keyString, Seq(JwtAlgorithm.HS512))
          case _ => throw new IllegalArgumentException("Only algorithms RS256,RS384,RS512,HS256,HS384,HS512 are supported.")
        }
      )
    }
  }

  /**
    * First get the claim without verification, to read the app name from the "iss" claim
    * and read the corresponding public key, then validate the token.
    * Return a Success only if all these steps are successful.
    */
  def validate(jwt: String, db: Database): Try[User] = {
    JwtJson.decodeJson(jwt, JwtOptions(signature = false)) flatMap { claim =>
      for {
        iss <- Try(issFromClaim(claim))
        username <- Try(userFromClaim(claim))
        (appId, keyString, algorithm) <- Try(appFromDatabase(iss))
        user <- Try(userFromDatabase(username, appId))
        _ <- validateToken(jwt, keyString, algorithm)
      } yield {
        user
      }
    }
  }

}
