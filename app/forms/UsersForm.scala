package forms

import models.User
import play.api.libs.json.JsValue



object UsersForm {
  /**
    * If appId is found in JSON, use it, otherwise can overwrite it with the *appId* parameter (defaults to 0).
    */
  def fromJson(j: JsValue, appId: Int = 0): User = {
    User(
      (j \ "appId").asOpt[Int] getOrElse appId,
      (j \ "username").as[String],
      (j \ "group").asOpt[String],
      (j \ "isActive").asOpt[Boolean] getOrElse false,
      (j \ "isAdmin").asOpt[Boolean] getOrElse false
    )
  }
}
