package auth

import play.api.libs.concurrent.Execution.Implicits._
import play.api.mvc.Results._
import play.api.mvc._
import play.api.http.HeaderNames._

import scala.concurrent.Future
import models.User


// N.B. "Unauthorized" (401) means "unauthenticated". "Forbidden" (403) means "not enought rights" (inferior role).

class AuthenticatedRequest[A](val user: User, request: Request[A]) extends WrappedRequest[A](request)


/**
  * Checks that the request's session JWT contains a known User.
  * If there is a User but he is not member of the facility, return status is Forbidden (403).
  * If there is no User, return status is Unauthorized (401).
  */
object AuthenticatedAction extends ActionBuilder[AuthenticatedRequest] {
  def getJWT[A](request: Request[A]) = {
    request.headers.get(AUTHORIZATION) match {
      case None => Future.successful(Unauthorized("No Authorization header"))
      case Some(jwt: String) => ???
    }
  }

  def invokeBlock[A](request: Request[A], block: AuthenticatedRequest[A] => Future[Result]) = {
//    request.jwtSession.getAs[JWTSession]("user") match {
//      case Some(user: JWTSession) if user.isFacility => block(new AuthenticatedRequest(user, request)).map(_.refreshJwtSession(request))
//      case Some(_) => Future.successful(Forbidden.refreshJwtSession(request))
//      case _ => Future.successful(Unauthorized)
//    }
    Future.successful(Ok)
  }
}