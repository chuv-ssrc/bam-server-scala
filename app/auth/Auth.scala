package auth

import java.io.File
import java.security.PublicKey

import pdi.jwt.{JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import models.User
import utils.Common._
import auth.RSA._


object Auth {

  def getPublicKey(appName: String, path: String = "resources/rsa_keys"): PublicKey = {
    findInTree(path, "cer") find {_.getName.contains(appName)} map {
      cert: File => readPublicKeyFromCertificate(cert.getPath)
    } getOrElse {
      findInTree(path, "pem") find {_.getName.contains(appName)} map {
        keyFile: File => readPublicKeyFromPemFile(keyFile.getPath)
      } getOrElse {
        throw new IllegalArgumentException("Could not find a suitable RSA public key in either .cer or .pem file.")
      }
    }
  }

  def validateToken(jwt: String, appName: String): Try[User] = {
    val publicKey = getPublicKey(appName)
    JwtJson.decodeJson(jwt, publicKey, Seq(JwtAlgorithm.RS256)) flatMap {
       claim => claim.asOpt[User] match {
         case None => Failure(new IllegalArgumentException("Could not cast JWT body/claim to case class User"))
         case Some(x) => Success(x)
       }
    }
  }

}
