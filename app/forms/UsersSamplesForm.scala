package forms

import models.UserSample
import play.api.libs.json.JsValue


object UsersSamplesForm {
  def fromJson(j: JsValue): UserSample = {
    UserSample(
      (j \ "username").as[String],
      (j \ "sample").as[String]
    )
  }
}
