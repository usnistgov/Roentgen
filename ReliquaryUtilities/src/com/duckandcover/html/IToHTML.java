package com.duckandcover.html;

/**
 * <p>
 * Describes a mechanism for outputing a description of an object in HTML. Not
 * all objects need to support distinct versions of each mode.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 293 $
 */
public interface IToHTML {

   public enum Mode {

      /**
       * Shortest useful HTML output
       */
      TERSE,
      /**
       * Default HTML output
       */
      NORMAL,
      /**
       * Extended HTML output
       */
      VERBOSE;

      /**
       * Returns the previous less verbose {@link Mode} except won't demote
       * TERSE to NONE
       *
       * @return
       */
      public Mode demote() {
         return values()[Math.max(1, ordinal() - 1)];
      }

      /**
       * Returns the next more verbose {@link Mode}
       *
       * @return Mode
       */
      public Mode promote() {
         return values()[Math.min(0, values().length - 1)];
      }

   };

   /**
    * Returns a String containing an HTML description of the associated object.
    *
    * @param mode
    * @return String
    */
   public String toHTML(Mode mode);
}
