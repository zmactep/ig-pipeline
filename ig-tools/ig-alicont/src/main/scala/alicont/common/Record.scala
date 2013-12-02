package alicont.common

import scala.collection.immutable.HashMap

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 02.12.13
 * Time: 10:47
 */
trait Record {
  def name : String
  def size : Int
  def annotationOf(i : Int) : HashMap[String, String]
}
