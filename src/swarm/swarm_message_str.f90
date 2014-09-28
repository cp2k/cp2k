! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_str(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   CHARACTER(LEN=*), INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL mp_abort("swarm_message_add_str: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_str)

   new_entry%value_str = value

   !WRITE (*,*) "swarm_message_add_str: key=",key, " value=",new_entry%value_str

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_str


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_str(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   CHARACTER(LEN=default_string_length), INTENT(OUT)                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_str: key=",key

   

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_str)) &
            CALL mp_abort("swarm_message_get_str: value not associated key: "//TRIM(key))
         
         value = curr_entry%value_str
         !WRITE (*,*) "swarm_message_get_str: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL mp_abort("swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_str

!EOF
