package models

import play.api.libs.json.Json


case class UserSample(
  username: String,    // The user identifier
  sample: String       // The sample identifier
)

object UserSample {
  implicit val userSampleFormat = Json.format[UserSample]
}
