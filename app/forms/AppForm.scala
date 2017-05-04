package forms

import models.App
import play.api.libs.json.JsValue


object AppForm {
  def fromJson(j: JsValue): App = {
    App(
      (j \ "iss").as[String],
      (j \ "keyFile").as[String],
      (j \ "description").asOpt[String]
    )
  }
}
