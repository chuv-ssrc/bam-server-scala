package forms

import models.User
import play.api.libs.json.JsValue



object UsersForm {
  def fromJson(j: JsValue, appId: Option[Int] = None): User = {
    User(
      (j \ "appId").asOpt[Int] orElse appId,
      (j \ "username").as[String],
      (j \ "group").asOpt[String],
      (j \ "isActive").asOpt[Boolean] getOrElse false,
      (j \ "isAdmin").asOpt[Boolean] getOrElse false
    )
  }
}
