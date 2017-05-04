package models

import play.api.libs.json._


case class User(
  appId: Option[Int],  // The client app identifier
  username: String,    // The user identifier, expecting a custom "name" claim in JWT (could use "sub" instead?)
  group: Option[String] = None,
  isActive: Boolean = false,
  isAdmin: Boolean = false
)

object User {
  implicit val userFormat = Json.format[User]
}
