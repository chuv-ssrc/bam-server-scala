package auth

import pdi.jwt.{JwtJson, JwtAlgorithm}
import scala.util.{Try, Success, Failure}
import play.api.libs.json._
import models.User



object Auth {

  def decodeJWT(jwt: String, secret: String): Try[JsObject] = {
    JwtJson.decodeJson(jwt, secret, Seq(JwtAlgorithm.RS256))
  }

}
