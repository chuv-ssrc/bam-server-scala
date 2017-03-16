package scala.auth

import pdi.jwt.{JwtJson, JwtAlgorithm}
import scala.util.{Try, Success, Failure}
import play.api.libs.json._

import org.scalatestplus.play._
import auth.Auth._



class AuthSpec extends PlaySpec {

  "`validateJWT`" should {

    "return the claim when the JWT is valid" in {
      decodeJWT("what?", "secret?") mustBe a[Success[_]]
    }

    "fail when the JWT is not valid" in {
      decodeJWT("what?", "secret?") mustBe a[Failure[_]]
    }

  }

}
