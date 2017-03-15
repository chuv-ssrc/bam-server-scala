package controllers

import javax.inject._
import play.api._
import play.api.mvc._


@Singleton
class HomeController @Inject() extends Controller {

 /**
   * Prints a message if the service works.
   */
  def index = Action {
    Ok("BAM server operational.")
  }

}
