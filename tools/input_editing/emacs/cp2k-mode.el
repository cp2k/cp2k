;;;;; Emacs major mode for cp2k input, written by Lianheng Tong
;;;;; Copyright (c) Lianheng Tong
;;;;; Last modify date: Saturday, 2014/01/25

;;;; Syntax highlighting of keywords
(defconst cp2k-font-lock-keywords
  (list
   ;; blocks
   '("^[ \t]*\\(&\\(END[ \t]*\\)?\\(\\sw+\\)\\)"
     (1 font-lock-function-name-face))
   ;; keywords
   '("^[ \t]*\\<\\([^&@$0-9]\\sw+\\)[ \t]+\\(\\sw+\\)?"
     (1 font-lock-keyword-face))
   ;; preprocessor type 1
   '("^[ \t]*\\(@\\(\\(I\\(F\\|NCLUDE\\)\\)\\|ENDIF\\)\\)"
     (1 font-lock-preprocessor-face))
   ;; preprocessor type 2
   '("^[ \t]*\\(@SET\\)[ \t]+\\(\\sw+\\)"
     (1 font-lock-preprocessor-face) (2 font-lock-variable-name-face nil t))
   ;; variables
   '("\\${?\\(\\sw+\\)}?"
     (1 font-lock-variable-name-face)))
  "Highlighting keywords for `cp2k-mode'.")

;;;; Syntax table
(defconst cp2k-mode-syntax-table
  (let ((st (make-syntax-table)))
    (modify-syntax-entry ?# "<" st)   ; #  is beg. comment
    (modify-syntax-entry ?\n ">" st)  ; \n is end. comment
    (modify-syntax-entry ?_ "w" st)   ; underscore is part of names
    (modify-syntax-entry ?\' "\"" st) ; string quote
    (modify-syntax-entry ?\" "\"" st) ; string quote
    (modify-syntax-entry ?\r "-" st)  ; return is whitespace
    (modify-syntax-entry ?+ "." st)   ; + is puntuation
    (modify-syntax-entry ?- "." st)   ; - is puntuation
    (modify-syntax-entry ?* "." st)   ; * is puntuation
    (modify-syntax-entry ?/ "." st)   ; / is puntuation
    (modify-syntax-entry ?= "." st)   ; = is puntuation
    (modify-syntax-entry ?\\ "\\" st) ; \ is escape char
    st)
  "Syntax table for `cp2k-mode'.")

;;;; Define keymap for the major mode
(defvar cp2k-mode-map
  (let ((map (make-sparse-keymap)))
    ;; Define mode specific key-bindings here
    (define-key map (kbd "C-j")     'newline-and-indent)
    (define-key map (kbd "C-c ;")   'comment-region)
    (define-key map (kbd "TAB")     'cp2k-indent-line)
    (define-key map (kbd "C-M-a")   'cp2k-beginning-of-block)
    (define-key map (kbd "C-M-e")   'cp2k-end-of-block)
    (define-key map (kbd "C-c C-c") 'outline-toggle-children)
    (define-key map (kbd "C-c C-a") 'show-all)
    (define-key map (kbd "C-c C-t") 'show-subtree)
    map)
  "Keymap for `cp2k-mode'.")

;;;; Syntax indentation
(defvar cp2k-indent
  2
  "standard indentation for in `cp2k-mode'")

(defconst cp2k-emptyline
  "^\\s-*$"
  "regexp matching an empty line in `cp2k-mode'")

(defconst cp2k-opening
  ;; note that this definition also matches cp2k-closing, so for
  ;; positive identification of cp2k-opening, one must use AND
  ;; construct with (not lookingat cp2k-closing).  One way to avoid
  ;; this is to construct search for line that does NOT contain END,
  ;; however, practical implementation of that leads to emacs
  ;; occasionally throwing the exception "regex stack overflow"...
  ;; The emacs regex search is pretty crudely implemented, and it is
  ;; far more robust to do a positive search, than negative
  ;; searches. So I opted for this more crude method here.
  "^[ \t]*\\(\\(&\\sw+\\)\\|\\(@IF\\>\\)\\)"
  "regexp matching lines starting (or closing --- due to
  limitations in regex implementation) a cp2k block in
  `cp2k-mode'")

(defconst cp2k-closing
  "^[ \t]*\\(\\(&END[ \t]*.*$\\)\\|\\(@ENDIF\\>\\)\\).*$"
  "regexp matching lines closing a cp2k block in `cp2k-mode'")

(defun cp2k-beginning-of-block (&optional not-set-mark)
  "move the cursor to the beginning of the block in `cp2k-mode'"
  (interactive)
  (if (not not-set-mark)
      (push-mark))
  (beginning-of-line)
  (let ((find_opening t)
        (in_nested nil)
        (initial_opening nil))
    ;; if looking at opening now, search for the parent opening
    (if (and (looking-at cp2k-opening)
             (not (looking-at cp2k-closing))
             (not (bobp)))
        (setq initial_opening t))
    (while (and find_opening (not (bobp)))
      (if (and (looking-at cp2k-opening)
               (not (looking-at cp2k-closing))
               (not in_nested)
               (not initial_opening))
          (setq find_opening nil)
        (progn
          (when (and (looking-at cp2k-opening)
                     (not (looking-at cp2k-closing)))
            (setq in_nested nil)
            (setq initial_opening nil))
          (forward-line -1)
          ;; if the line is a closing statement then
          ;; skip to the line of the opening statement
          (if (looking-at cp2k-closing)
              (progn
                (setq in_nested t)
                (cp2k-beginning-of-block t))))))))

(defun cp2k-end-of-block (&optional not-set-mark)
  "move the cursor to the ending of the block in `cp2k-mode'"
  (interactive)
  (if (not not-set-mark)
      (push-mark))
  (beginning-of-line)
  (let ((find_closing t)
        (in_nested nil)
        (initial_closing nil))
    ;; if looking at closing now, search for the parent closing
    (if (and (looking-at cp2k-closing)
             (not (eobp)))
        (setq initial_closing t))
    (while (and find_closing (not (eobp)))
      (if (and (looking-at cp2k-closing)
               (not in_nested)
               (not initial_closing))
          (setq find_closing nil)
        (progn
          (when (looking-at cp2k-closing)
            (setq in_nested nil)
            (setq initial_closing nil))
          (forward-line 1)
          ;; if the line is a opening statement then
          ;; skip to the line of the closing statement
          (if (and (looking-at cp2k-opening)
                   (not (looking-at cp2k-closing)))
              (progn
                (setq in_nested t)
                (cp2k-end-of-block t))))))))

(defun cp2k-forward-one-line ()
  "move the cursor foward one line, ignore empty or comment lines, in `cp2k-mode'"
  (beginning-of-line)
  (forward-line 1)
  (while (and (or (looking-at cp2k-emptyline)
                  (looking-at "^[ \t]*#"))
              (not (eobp)))
    (forward-line 1)))

(defun cp2k-backward-one-line ()
  "move the cursor foward one line, ignore empty or comment lines, in `cp2k-mode'"
  (beginning-of-line)
  (forward-line -1)
  (while (and (or (looking-at cp2k-emptyline)
                  (looking-at "^[ \t]*#"))
              (not (bobp)))
    (forward-line -1)))

(defun cp2k-left-of-point-is-empty ()
  "check if the left of the point is only white space in `cp2k-mode'"
  (let (not-empty)
    (save-excursion
      (setq not-empty
            (re-search-backward "[^ \t]" (line-beginning-position) t 1)))
    (not not-empty)))

(defun cp2k-indent-line ()
  "indent current line in `cp2k-mode'"
  (interactive)
  (let ((indent 0)
        (position (point-marker)))
    ;; record initial indentation position
    (back-to-indentation)
    (indent-line-to
     (max
      0
      (catch 'indentation
        (save-excursion
          (beginning-of-line)
          (if (bobp) (throw 'indentation 0))
          ;; get the block indentation
          (save-excursion
            (cp2k-beginning-of-block t)
            (setq indent (current-indentation)))
          ;; look at the current line
          (if (looking-at cp2k-closing)
              (throw 'indentation indent))
          ;; move up one non-empty line
          (cp2k-backward-one-line)
          (setq indent (current-indentation))
          (if (and (looking-at cp2k-opening)
                   (not (looking-at cp2k-closing)))
              (setq indent (+ indent cp2k-indent)))
          (throw 'indentation indent)))))
    ;; move cursor to indentation point if in the beginning white
    ;; space otherwise leave unchanged
    (goto-char position)
    (if (cp2k-left-of-point-is-empty)
        (back-to-indentation))
    (set-marker position nil)))

;;;; Define mode hook
(defvar cp2k-mode-hook nil)

;;;; Entry function
(define-derived-mode cp2k-mode fundamental-mode "cp2k"
  "Major mode for editing cp2k input. Copyright (c) Lianheng Tong 2013/12/06"
  :syntax-table cp2k-mode-syntax-table
  ;; setq-local is undefined in emacs versions prior 24.3, it is
  ;; equivalent to make-local-variable, and followed by setq, and this
  ;; is what has been used here. This works for older versions of
  ;; emacs.
  (make-local-variable 'comment-start)
  (setq comment-start "# ")
  (make-local-variable 'comment-start-skip)
  (setq comment-start-skip "#+\\s-*")
  ;; see emacs manual on font-lock-defaults, nil here means also
  ;; font-lock comments and strings, and t means the value of
  ;; font-lock-keywords-case-fold-search is set to non-nil, allowing
  ;; case insensitive search
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(cp2k-font-lock-keywords nil t))
  ;; whether the indentation key-word search case-sensitive or not (in
  ;; my implementation, using looking-at function) is effected only by
  ;; the variable case-fold-search, and is uneffected by
  ;; font-lock-keywords-case-fold-search. So turn case-fold-search on
  ;; just in case
  (make-local-variable 'case-fold-search)
  (setq case-fold-search t)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'cp2k-indent-line))

;;;; Add to autoload alist
(add-to-list 'auto-mode-alist '("\\.cp2kin\\'" . cp2k-mode))

;;;; Setup outline mode
(when (require 'outline nil 'noerror)
  (add-hook 'cp2k-mode-hook 'outline-minor-mode)
  ;; in outline mode, the level of header depends on the length of
  ;; match, which suits our purposes quite well
  (setq outline-regexp "[ \t]*\\(&\\|@\\(IF\\|EN\\)\\)"))

;;;; Last line
(provide 'cp2k-mode)
