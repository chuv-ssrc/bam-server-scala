package models

import play.api.libs.json._


case class App(
  iss: String,  // as in JWT "iss" claim
  key: String,  // path to the RSA public key file
  description: Option[String] = None
)

object App {
  implicit val appFormat = Json.format[App]
}

