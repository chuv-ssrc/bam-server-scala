package models

import play.api.libs.json._


case class User(
  username: String,
  role: String,
  samples: List[String]
) {
  implicit val userFormat = Json.format[User]
}
