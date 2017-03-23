package auth

import pdi.jwt.{JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import models.User


object Auth {

//  def decodeJWT_RSA(jwt: String, secret: String): Try[JsObject] = {
//    JwtJson.decodeJson(jwt, secret, Seq(JwtAlgorithm.RS256))
//  }

}
