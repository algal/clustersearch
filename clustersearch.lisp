(defun pdf-valid-p (pdf)
  (equal 1 (reduce #'+ pdf)))

(defun cdf<-pdf (pdf)
  (loop for p in pdf
        with running-sum = 0
        collecting (incf running-sum p)))

(defun random-real ()
  "Returns real in [0,1)

   See: http://www.cs.cmu.edu/Groups/AI/html/cltl/clm/node133.html"
  (let ((maxint 1000))
    (float (/ (random maxint) maxint))))

(defun random-with-cdf (cdf)
  "Returns random non-negative integer.

   Chosen based on specified cumulative distribution."
  (let ((r (random-real)))
    (loop for bucket-max in cdf
          for chosen from 0
          thereis (when (< r bucket-max) chosen))))

(defun random-with-pdf (pdf)
  (random-with-cdf (cdf<-pdf pdf)))
