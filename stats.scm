; minmax comp is comparator function (< > <= >=)
(define (minmax ls comp)
  (define n (length ls))
  (define (iter ls k m)
    (if (> k n)
      m
      (iter (cdr ls) (+ k 1) (if (comp (car ls) m) (car ls) m))))
  (iter ls 1 (car ls)))

; mean average
(define (mean ls)
  (let ((n (* 1.0 (length ls))))
    (/ (apply + ls) n)))

; return nth element of ls - necessary auxiliary to median procedure
(define (nth ls n)
  (define (iter ls k)
    (if (= n k)
      (car ls)
      (iter (cdr ls) (+ k 1))))
  (iter ls 1))

; median average
(define (median ls)
  (let ((sorted-ls (sort ls <)) (m (ceiling (/ (length ls) 2))))
    (if (even? (length ls))
      (/ (+ (nth sorted-ls m) (nth sorted-ls (+ m 1))) 2.0)
      (nth sorted-ls m))))

; construct list of deviations about a central tendency
(define (dev ls m) (map (lambda (x) (- x m)) ls))

; list of absolute values of deviations about central tendency
(define (adev ls m) (map abs (dev ls m)))

; compute population variance of list
(define (population-variance ls)
  (let ((n (length ls)))
    (/ (apply + (map square (adev ls (mean ls)))) n)))

; population standard deviation
(define (population-standard-deviation ls) (sqrt (population-variance ls)))

; compute sample variance of list
(define (sample-variance ls)
  (let ((n (length ls)))
    (/ (apply + (map square (adev ls (mean ls)))) (- n 1))))

; sample standard deviation
(define (sample-standard-deviation ls) (sqrt (sample-variance ls)))

; median absolute deviation
(define (median-deviation ls) (median (adev ls (median ls))))

; mean absolute deviation
(define (mean-deviation ls) (mean (adev ls (mean ls))))

; cube auxiliary for skewness procedures
(define (cube x) (* x x x))

; population skewness
(define (population-skewness ls)
  (let ((xdev (dev ls (mean ls))) 
        (n (length ls))
        (std population-standard-deviation))
    (/ (/ (apply + (map cube xdev)) n) (cube (std ls)))))

; sample skewness
(define (sample-skewness ls)
  (let ((xdev (dev ls (mean ls))) 
        (n (- (length ls) 1))
        (std sample-standard-deviation))
    (/ (/ (apply + (map cube xdev)) n) (cube (std ls)))))

; non-parametric skew
(define (non-parametric-skew ls)
  (/ (- (mean ls) (median ls)) (population-standard-deviation ls)))

; covariance of two data sets
(define (covariance lsx lsy)
  (let ((xdev (dev lsx (mean lsx))) 
        (ydev (dev lsy (mean lsy)))
        (n (min (length lsx) (length lsy))))
    (define (iter xdev ydev s)
      (if (or (null? xdev) (null? ydev))
        (/ s n)
        (iter (cdr xdev) (cdr ydev) (+ s (* (car xdev) (car ydev))))))
    (iter xdev ydev 0.0)))

; correlation coefficient of two data sets
(define (correlation-coefficient lsx lsy)
  (let ((std population-standard-deviation))
    (/ (covariance lsx lsy) (* (std lsx) (std lsy)))))

; data set generation procedures below

; uniform random distribution of float
(define (uniform n)
  (define (iter ls k)
    (if (> k n)
      ls
      (iter (cons (random 1.0) ls) (+ k 1))))
  (iter () 1))

; normal distribution of float
(define (normal n)
  (define pi (* 4.0 (atan 1.0)))
  (define mu (/ n 2.0))
  (define sd (/ n 6.024))
  (define (iter ls k)
    (if (> k n)
      ls
      (let ((r 
              (*
              (/ 1. (sqrt (* 2. pi (square sd))))
              (exp (/ (* -1. (square (- k mu)))
                      (* 2. (square sd)))))))
      (iter (cons (* sd r) ls) (+ k 1)))))
  (iter '() 1))

(define (normal-random n)
  (define pi (* 4.0 (atan 1.0)))
  (define mu (/ n 2.0))
  (define sd (/ n 6.024))
  (define (iter ls k)
    (if (> k n)
      ls
      (let ((r 
              (*
              (/ 1. (sqrt (* 2. pi (square sd))))
              (exp (/ (* -1. (square (- (random n) mu)))
                      (* 2. (square sd)))))))
      (iter (cons (* sd r) ls) (+ k 1)))))
  (iter '() 1))

(define (box-muller n)
  (define pi (* 4.0 (atan 1.0)))
  (define (iter ls k)
    (if (> k n)
      ls
      (if (odd? k)
        (let ((z
              (* 
                (sqrt (* -2.0 (log (random 1.0))))
                (cos (* 2.0 pi (random 1.0))))))
        (iter (cons z ls) (+ k 1)))
        (let ((z
              (* 
                (sqrt (* -2.0 (log (random 1.0))))
                (sin (* 2.0 pi (random 1.0))))))
        (iter (cons z ls) (+ k 1))))))
  (iter () 1))

(define (randnorm n)
  (define ls (box-muller n))
  (define a (minmax ls <))
  (define b (minmax ls >))
  (define range (- b a))
  (define c
    (if (> (* -1 a) b)
      a
      b))
  (define r (* 2 c))
  (define d (/ (- r range) 2.0))
  (define (iter ls sls)
    (if (null? ls)
      sls
      (let ((x (/ (+ (- (car ls) a) d) r)))
        (iter (cdr ls) (cons x sls)))))
  (iter ls ()))

(define (normal-int n)
  (define ls (randnorm n))
  (define (iter ls-in ls-out)
    (if (null? ls-in)
      ls-out
      (let ((x (inexact->exact (round (* 100 (car ls-in))))))
      (iter (cdr ls-in) (cons x ls-out)))))
  (iter ls ()))
