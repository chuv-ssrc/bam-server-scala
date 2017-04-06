package models

import play.api.libs.json._


case class User(
  name: String, // The user identifier, expecting a custom "name" claim in JWT (could use "sub" instead?)
  appId: Option[Int],  // The client app identifier
  isActive: Boolean = false,
  isAdmin: Boolean = false
)

object User {
  implicit val userFormat = Json.format[User]
}